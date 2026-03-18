#!/usr/bin/env python3
"""
graph_hgt_pipeline.py

END-TO-END ANALYSIS PIPELINE (STUBS) for the *already-pruned undirected* cross-species protein graph.

Input format (CSV/TSV), columns:
  u, v, shared_kmers, jaccard, species_u, species_v

Example row:
  GCF_000005845.2_ASM584v2|NP_414543.1,
  GCF_000240185.1_ASM24018v2|YP_005225012.1,
  544,0.505107,Escherichia coli,Klebsiella pneumoniae

Assumptions:
  - Edges are already undirected and pruned (degree explosion + per species-pair Jaccard percentile).
  - Edges are cross-species (species_u != species_v).
  - File size ~10–20MB; feasible to load fully into memory.

Outputs (recommended deliverables):
  1) edge_features.tsv
  2) protein_features.tsv
  3) component_features.tsv
  4) hgt_candidates.tsv
  5) (optional) high_z_subgraph_edges.tsv for "bridge gene" analyses

Key "additions" included:
  - Robust edge surprise (z_robust) via median + MAD per species-pair.
  - Component-level HGT signal: species entropy + high-z concentration.

Implementation notes:
  - Use networkx for graph algorithms (components, clustering coeff, betweenness).
  - For speed, compute betweenness on a reduced high-z subgraph or per-component.

Related: upstream candidate generator script you used earlier:
  /mnt/data/kmer_candidates_from_faa.py  :contentReference[oaicite:0]{index=0}
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

# Optional: make networkx a dependency for graph algorithms
try:
    import networkx as nx
except ImportError:  # pragma: no cover
    nx = None


import time
from contextlib import contextmanager

@contextmanager
def stage(name: str):
    """Context manager to time a pipeline stage."""
    t0 = time.perf_counter()
    print(f"[STAGE] {name} ...")
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        print(f"[DONE ] {name} in {dt:.3f}s")


def log_counts(edges=None, G=None, pair_stats=None, prefix: str = "[INFO]"):
    """Optional: print quick counts for sanity."""
    parts = []
    if edges is not None:
        parts.append(f"edges={len(edges)}")
    if G is not None:
        parts.append(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    if pair_stats is not None:
        parts.append(f"species_pairs={len(pair_stats)}")
    if parts:
        print(f"{prefix} " + " | ".join(parts))

# =========================
# Data structures
# =========================

@dataclass(frozen=True)
class Edge:
    """Undirected edge record."""
    u: str
    v: str
    shared_kmers: int
    jaccard: float
    species_u: str
    species_v: str

    def species_pair(self) -> Tuple[str, str]:
        """Return unordered species-pair key (A,B) with A <= B."""
        a, b = self.species_u, self.species_v
        return (a, b) if a <= b else (b, a)

    def node_pair(self) -> Tuple[str, str]:
        """Return unordered node-pair key (u,v) with u <= v."""
        a, b = self.u, self.v
        return (a, b) if a <= b else (b, a)


@dataclass
class PairRobustStats:
    """Robust baseline stats for a species-pair."""
    n: int
    median: float
    mad: float
    scale: float  # 1.4826*MAD + eps


@dataclass
class EdgeFeatures:
    """Computed edge-level features."""
    u: str
    v: str
    species_u: str
    species_v: str
    shared_kmers: int
    jaccard: float
    z_robust: Optional[float]  # None if insufficient data
    pair_median: Optional[float]
    pair_mad: Optional[float]
    pair_n: int


@dataclass
class NodeFeatures:
    """Computed node-level features."""
    u: str
    species: str
    component_id: int

    # Connectivity
    deg_xsp: int
    n_species: int
    top1_jac: float
    top5_mean_jac: float
    top1_shared: int
    top5_mean_shared: float

    # Surprise
    max_z: Optional[float]
    top5_mean_z: Optional[float]

    # Specificity / mixing
    max_species_fraction: float
    species_entropy: float
    species_entropy_norm: float
    n_eff_species: float
    participation_coeff: float

    # Bridge-ish topology
    clustering_coeff: Optional[float]
    betweenness: Optional[float]


@dataclass
class ComponentFeatures:
    """Computed component-level features."""
    component_id: int
    size: int
    num_edges: int
    n_species: int
    species_entropy: float
    species_entropy_norm: float
    high_z_frac: Optional[float]  # None if z unavailable
    max_z: Optional[float]
    top10_mean_z: Optional[float]


# =========================
# I/O
# =========================

def sniff_delimiter(path: Path) -> str:
    """
    Try to detect delimiter (',' or '\\t') from the first non-empty line.
    """
    with open(path, "r", encoding="utf-8", newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "\t" in line and "," not in line:
                return "\t"
            if "," in line and "\t" not in line:
                return ","
            # If both exist, fallback to CSV logic:
            return ","
    return ","


def read_edges(path: Path) -> List[Edge]:
    """
    Read edge list from CSV/TSV with header:
      u, v, shared_kmers, jaccard, species_u, species_v

    Returns a list of Edge records.

    NOTE: This assumes the file represents undirected edges already.
    If duplicates exist (u,v) and (v,u), you should dedupe (see dedupe_edges()).
    """
    delim = sniff_delimiter(path)
    edges: List[Edge] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter=delim)
        required = {"u", "v", "shared_kmers", "jaccard", "species_u", "species_v"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}. Found: {r.fieldnames}")

        for row in r:
            edges.append(
                Edge(
                    u=row["u"],
                    v=row["v"],
                    shared_kmers=int(row["shared_kmers"]),
                    jaccard=float(row["jaccard"]),
                    species_u=row["species_u"],
                    species_v=row["species_v"],
                )
            )
    return edges


def dedupe_edges(edges: List[Edge]) -> List[Edge]:
    """
    Dedupe duplicates for an undirected graph.

    Policy:
      For each unordered protein pair {u, v},
      keep exactly one edge — the one with the maximum Jaccard similarity.
    """
    best: Dict[Tuple[str, str], Edge] = {}
    for e in edges:
        key = e.node_pair()
        prev = best.get(key)
        if prev is None:
            best[key] = e
        else:
            if e.jaccard > prev.jaccard:
                best[key] = e
    return list(best.values())


# =========================
# Robust z-score per species pair (addition #1)
# =========================

def median(xs: Sequence[float]) -> float:
    """
    Compute the median of a non-empty sequence of numbers.

    Policy:
      - Sort once (O(n log n)), which is fine given per-species-pair sizes.
      - If n is odd: middle element
      - If n is even: average of the two middle elements
    """
    if not xs:
        raise ValueError("median() arg is an empty sequence")

    xs_sorted = sorted(xs)
    n = len(xs_sorted)
    mid = n // 2

    if n % 2 == 1:
        return xs_sorted[mid]
    else:
        return 0.5 * (xs_sorted[mid - 1] + xs_sorted[mid])


def mad(xs: Sequence[float], med: float) -> float:
    """
    Compute the Median Absolute Deviation (MAD) from a given median.

    MAD = median(|x - med|)

    Notes:
      - Uses absolute deviations from the *already computed* median.
      - Does not rescale by 1.4826 here (that belongs to the caller).
    """
    if not xs:
        raise ValueError("mad() arg is an empty sequence")

    abs_devs = [abs(x - med) for x in xs]
    return median(abs_devs)



def compute_pair_robust_stats(
    edges: List["Edge"],
    weight: str = "jaccard",
    min_pair_edges: int = 30,
    eps: float = 1e-9,
    shared_transform: str = "log1p",  # "raw" or "log1p"
) -> Dict[Tuple[str, str], "PairRobustStats"]:
    """
    For each unordered species pair (A,B), compute robust baseline stats over edge weights.

    Robust baseline:
      m = median(w)
      MAD = median(|w - m|)
      scale = 1.4826*MAD + eps

    Parameters:
      weight:
        - 'jaccard' (recommended)
        - 'shared'  (optionally transformed; see shared_transform)
      min_pair_edges:
        - Used later when deciding whether to trust z-scores.
        - Here we still compute stats for all pairs.
      eps:
        - Added to scale so z-scores don't divide by zero when MAD=0.
      shared_transform:
        - 'log1p' makes shared_kmers less dominated by very large proteins.
        - 'raw' uses shared_kmers as-is.

    Returns:
      dict: (A,B) -> PairRobustStats
    """
    by_pair: Dict[Tuple[str, str], List[float]] = defaultdict(list)

    if weight not in {"jaccard", "shared"}:
        raise ValueError("weight must be 'jaccard' or 'shared'")
    if shared_transform not in {"raw", "log1p"}:
        raise ValueError("shared_transform must be 'raw' or 'log1p'")

    # 1) bucket weights by unordered species pair
    if weight == "jaccard":
        for e in edges:
            by_pair[e.species_pair()].append(e.jaccard)
    else:
        if shared_transform == "log1p":
            for e in edges:
                by_pair[e.species_pair()].append(math.log1p(e.shared_kmers))
        else:
            for e in edges:
                by_pair[e.species_pair()].append(float(e.shared_kmers))

    # 2) compute robust stats per bucket
    stats: Dict[Tuple[str, str], PairRobustStats] = {}
    for key, ws in by_pair.items():
        m = median(ws)
        d = mad(ws, m)
        s = 1.4826 * d + eps
        stats[key] = PairRobustStats(n=len(ws), median=m, mad=d, scale=s)

    return stats


def print_pair_stats_sanity_table(
    pair_stats: Dict[Tuple[str, str], "PairRobustStats"],
    min_pair_edges_for_z: int = 30,
    top_n: int = 15,
) -> None:
    """
    Print a quick sanity table to visually verify species-pair baselines.

    Shows the top-N species pairs by edge count, with:
      - n edges
      - median weight
      - MAD
      - scale (1.4826*MAD)
      - "ok_for_z" flag (n >= min_pair_edges_for_z)

    Call this right after compute_pair_robust_stats().
    """
    rows = []
    for (a, b), st in pair_stats.items():
        ok = "YES" if st.n >= min_pair_edges_for_z else "NO"
        rows.append((st.n, a, b, st.median, st.mad, st.scale, ok))

    rows.sort(reverse=True, key=lambda t: t[0])
    rows = rows[:top_n]

    print("\n[PAIR-STATS] Robust baselines per species-pair (top by edge count)")
    print(f"{'n':>7}  {'species_A':<24}  {'species_B':<24}  {'median':>10}  {'MAD':>10}  {'scale':>10}  {'ok_for_z':>8}")
    print("-" * (7 + 2 + 24 + 2 + 24 + 2 + 10 + 2 + 10 + 2 + 10 + 2 + 8))
    for n, a, b, med, d, s, ok in rows:
        print(f"{n:7d}  {a:<24.24}  {b:<24.24}  {med:10.6f}  {d:10.6f}  {s:10.6f}  {ok:>8}")
    print()


def compute_edge_features(
    edges: List[Edge],
    pair_stats: Dict[Tuple[str, str], PairRobustStats],
    weight: str = "jaccard",
    min_pair_edges_for_z: int = 30,
) -> List[EdgeFeatures]:
    """
    Compute edge-level features, including robust z-score per species-pair.

    z(e) = (w(e) - median_{A,B}) / scale_{A,B}

    If a species-pair has fewer than min_pair_edges_for_z edges,
    z_robust is set to None (unreliable baseline).

    Parameters:
      edges: list of undirected Edge objects
      pair_stats: output of compute_pair_robust_stats()
      weight: 'jaccard' or 'shared'
      min_pair_edges_for_z: minimum edges needed to trust z-score

    Returns:
      List[EdgeFeatures]
    """
    out: List[EdgeFeatures] = []

    if weight not in {"jaccard", "shared"}:
        raise ValueError("weight must be 'jaccard' or 'shared'")

    for e in edges:
        key = e.species_pair()
        st = pair_stats.get(key)

        # choose edge weight
        if weight == "jaccard":
            w = e.jaccard
        else:
            w = float(e.shared_kmers)  # (log transform already handled upstream if desired)

        if st is None:
            # should not happen, but keep code robust
            z = None
            pair_median = None
            pair_mad = None
            pair_n = 0
        else:
            pair_median = st.median
            pair_mad = st.mad
            pair_n = st.n

            if st.n < min_pair_edges_for_z:
                z = None
            else:
                # robust z-score
                z = (w - st.median) / st.scale

        out.append(
            EdgeFeatures(
                u=e.u,
                v=e.v,
                species_u=e.species_u,
                species_v=e.species_v,
                shared_kmers=e.shared_kmers,
                jaccard=e.jaccard,
                z_robust=z,
                pair_median=pair_median,
                pair_mad=pair_mad,
                pair_n=pair_n,
            )
        )

    return out



def write_edge_features(path: Path, feats: List[EdgeFeatures]) -> None:
    """
    Write edge_features.tsv with robust z and per-pair baseline fields.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "u", "v", "species_u", "species_v",
            "shared_kmers", "jaccard",
            "z_robust", "pair_median", "pair_mad", "pair_n"
        ])
        for e in feats:
            w.writerow([
                e.u, e.v, e.species_u, e.species_v,
                e.shared_kmers, f"{e.jaccard:.6f}",
                "" if e.z_robust is None else f"{e.z_robust:.6f}",
                "" if e.pair_median is None else f"{e.pair_median:.6f}",
                "" if e.pair_mad is None else f"{e.pair_mad:.6f}",
                e.pair_n,
            ])


# =========================
# Graph construction
# =========================

def build_graph(edges: List[Edge], use_networkx: bool = True):
    """
    Construct the undirected graph object.

    If use_networkx=True, returns a networkx.Graph with:
      node attrs: species
      edge attrs: shared_kmers, jaccard
    """
    if use_networkx and nx is None:
        raise RuntimeError("networkx not installed; set use_networkx=False or install it.")

    if use_networkx:
        G = nx.Graph()
        for e in edges:
            # Node attrs
            if e.u not in G:
                G.add_node(e.u, species=e.species_u)
            if e.v not in G:
                G.add_node(e.v, species=e.species_v)
            # Edge attrs
            G.add_edge(e.u, e.v, shared_kmers=e.shared_kmers, jaccard=e.jaccard)
        return G

    # TODO: implement your own adjacency list structure if you want no networkx
    #raise NotImplementedError


def attach_edge_z_to_graph(G, edge_feats: List[EdgeFeatures]) -> None:
    """
    Attach robust z-score to each graph edge.

    For every EdgeFeatures entry (u, v), set:
        G[u][v]["z_robust"] = z_robust

    Assumptions:
      - G is an undirected graph
      - G contains exactly the same undirected edges as edge_feats
      - edge_feats already reflect deduplicated edges

    Raises:
      KeyError if an edge in edge_feats does not exist in G
    """
    missing = 0

    for ef in edge_feats:
        u, v = ef.u, ef.v

        if not G.has_edge(u, v):
            # This should never happen if the pipeline is correct
            missing += 1
            continue

        # Attach z (may be None)
        G[u][v]["z_robust"] = ef.z_robust

        # Optional but useful: attach pair baseline info for debugging
        G[u][v]["pair_n"] = ef.pair_n
        G[u][v]["pair_median"] = ef.pair_median
        G[u][v]["pair_mad"] = ef.pair_mad

    if missing > 0:
        raise KeyError(
            f"attach_edge_z_to_graph: {missing} edges from edge_feats "
            f"were not found in the graph"
        )



# =========================
# Node features (connectivity, specificity, surprise, bridge-ness)
# =========================

from typing import Sequence, Dict, List, Optional
import math

def safe_mean(xs: Sequence[float]) -> float:
    """
    Mean with empty handling.

    Returns 0.0 if xs is empty (useful for features like top5_mean_z when no z exists).
    """
    if not xs:
        return 0.0
    return sum(xs) / len(xs)


def topk_mean(xs: Sequence[float], k: int) -> float:
    """
    Mean of the top-k values in xs.

    Efficient enough for per-node neighbor lists (degrees are usually capped by pruning).
    If xs has fewer than k elements, uses all of them.
    Returns 0.0 on empty.
    """
    if not xs:
        return 0.0
    if k <= 0:
        raise ValueError("k must be positive")
    # Sorting is fine since per-node degrees are limited by your pruning.
    ys = sorted(xs, reverse=True)
    ys = ys[:k]
    return sum(ys) / len(ys)


def entropy_from_counts(counts: Dict[str, int]) -> float:
    """
    Shannon entropy (natural log) over categorical counts.

    H = - Σ_s p_s log(p_s)
    where p_s = count_s / total.

    Returns 0.0 if total=0 or only one category.
    """
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c <= 0:
            continue
        p = c / total
        h -= p * math.log(p)
    return h


def participation_coefficient(counts: Dict[str, int], deg: int) -> float:
    """
    Participation coefficient treating each species as a 'community':
      P = 1 - Σ_s (c_s/deg)^2

    Interpretation:
      - P ≈ 0: neighbors concentrated in one species
      - P → 1: neighbors spread across many species

    Returns 0.0 if deg=0.
    """
    if deg <= 0:
        return 0.0
    s = 0.0
    for c in counts.values():
        frac = c / deg
        s += frac * frac
    return 1.0 - s


def compute_node_features(
    G,
    component_id_of: Dict[str, int],
    z_attr: str = "z_robust",
    compute_clustering: bool = True,
    compute_betweenness: bool = True,
    betweenness_mode: str = "high_z",  # 'full', 'per_component', 'high_z'
    z_threshold_for_high_subgraph: float = 3.0,
) -> List[NodeFeatures]:
    """
    Compute per-node features:

    Connectivity:
      deg_xsp, n_species, top1/top5 on jaccard/shared.

    Specificity:
      max_species_fraction
      species entropy (and normalized entropy)
      n_eff = exp(H)
      participation coefficient

    Surprise:
      max_z, top5_mean_z (if z available)

    Bridge-ness:
      clustering coeff (cheap)
      betweenness centrality (expensive; recommended to run on high-z subgraph or per-component)

    betweenness_mode:
      - 'full': betweenness on full graph (can be heavy)
      - 'per_component': compute betweenness per connected component
      - 'high_z': compute betweenness only on subgraph edges with z >= threshold
   Compute per-node features.

    Requirements in G:
      - node attr: 'species'
      - edge attrs: 'jaccard', 'shared_kmers', and optionally z_attr (may be None)

    Returns:
      List[NodeFeatures] (one per node)
    """
    if nx is None:
        raise RuntimeError("networkx is required for compute_node_features()")

    if betweenness_mode not in {"full", "per_component", "high_z"}:
        raise ValueError("betweenness_mode must be one of: full, per_component, high_z")

    # ---- 1) clustering coefficient (cheap) ----
    clustering: Optional[Dict[str, float]] = None
    if compute_clustering:
        # For undirected graphs, nx.clustering is straightforward.
        clustering = nx.clustering(G)

    # ---- 2) betweenness centrality (expensive) ----
    betweenness: Optional[Dict[str, float]] = None
    if compute_betweenness:
        if betweenness_mode == "full":
            betweenness = nx.betweenness_centrality(G, normalized=True)

        elif betweenness_mode == "per_component":
            betweenness = {}
            # compute per connected component to reduce cost
            for nodes in nx.connected_components(G):
                if len(nodes) <= 2:
                    for u in nodes:
                        betweenness[u] = 0.0
                    continue
                H = G.subgraph(nodes)
                b = nx.betweenness_centrality(H, normalized=True)
                betweenness.update(b)


        else:  # "high_z"

            print(f"[BETWEENNESS] mode=high_z z_attr={z_attr} z_threshold={z_threshold_for_high_subgraph}")

            # Subgraph with only edges that have z >= threshold

            H = nx.Graph()

            # Preserve node attrs only for nodes that end up in H (avoid adding all nodes up front).

            # This is a BIG speed/memory win.

            kept_edges = 0

            t0 = time.perf_counter()

            last_log = t0

            scanned = 0

            total_edges = G.number_of_edges()

            for u, v, data in G.edges(data=True):

                scanned += 1

                z = data.get(z_attr)

                if z is None or z < z_threshold_for_high_subgraph:

                    # periodic progress log while scanning

                    now = time.perf_counter()

                    if now - last_log >= 2.0:  # every ~2 seconds

                        pct = 100.0 * scanned / max(1, total_edges)

                        print(
                            f"[BETWEENNESS] scanning edges: {scanned}/{total_edges} ({pct:.1f}%) kept_edges={kept_edges}")

                        last_log = now

                    continue

                # add nodes lazily

                if u not in H:
                    H.add_node(u, **G.nodes[u])

                if v not in H:
                    H.add_node(v, **G.nodes[v])

                H.add_edge(u, v)

                kept_edges += 1

                now = time.perf_counter()

                if now - last_log >= 2.0:
                    pct = 100.0 * scanned / max(1, total_edges)

                    print(f"[BETWEENNESS] scanning edges: {scanned}/{total_edges} ({pct:.1f}%) kept_edges={kept_edges}")

                    last_log = now

            dt_scan = time.perf_counter() - t0

            print(
                f"[BETWEENNESS] high-z subgraph built in {dt_scan:.2f}s: H_nodes={H.number_of_nodes()} H_edges={H.number_of_edges()}")

            if kept_edges == 0:

                betweenness = {u: 0.0 for u in G.nodes()}

            else:

                # Log component size distribution in H (critical for understanding runtime)

                t1 = time.perf_counter()

                comps = [len(c) for c in nx.connected_components(H)]

                comps.sort(reverse=True)

                largest = comps[0] if comps else 0

                ncomp = len(comps)

                med = comps[len(comps) // 2] if comps else 0

                print(f"[BETWEENNESS] H components: n={ncomp} largest={largest} median={med}")

                # If H is still too big, approximate betweenness.

                # Rule of thumb: exact betweenness gets painful when largest component is big.

                APPROX_IF_LARGEST_GE = 800  # tune this

                APPROX_K = 200  # number of sampled sources (Brandes approx)

                if largest >= APPROX_IF_LARGEST_GE:

                    print(f"[BETWEENNESS] H largest component={largest} -> using APPROX betweenness (k={APPROX_K})")

                    b = nx.betweenness_centrality(H, k=APPROX_K, seed=42, normalized=True)

                else:

                    print(f"[BETWEENNESS] computing EXACT betweenness on H ...")

                    b = nx.betweenness_centrality(H, normalized=True)

                dt_bw = time.perf_counter() - t1

                print(f"[BETWEENNESS] betweenness done in {dt_bw:.2f}s")

                # Fill missing nodes with 0 on full graph node set

                betweenness = {u: b.get(u, 0.0) for u in G.nodes()}

    # ---- 3) per-node aggregation ----
    feats: List[NodeFeatures] = []
    for u in G.nodes():
        sp_u = G.nodes[u].get("species", "")

        # neighbors and edge attributes
        neigh_species = Counter()
        jac_list: List[float] = []
        shared_list: List[float] = []
        z_list: List[float] = []

        for v, edata in G[u].items():
            # neighbor species
            sp_v = G.nodes[v].get("species", "")
            neigh_species[sp_v] += 1

            # edge weights
            j = edata.get("jaccard")
            if j is not None:
                jac_list.append(float(j))

            sh = edata.get("shared_kmers")
            if sh is not None:
                shared_list.append(float(sh))

            z = edata.get(z_attr)
            if z is not None:
                z_list.append(float(z))

        deg = len(G[u])
        n_species = len(neigh_species)

        # top stats
        top1_jac = max(jac_list) if jac_list else 0.0
        top5_mean_jac = topk_mean(jac_list, 5)

        top1_shared = int(max(shared_list)) if shared_list else 0
        top5_mean_shared = topk_mean(shared_list, 5)

        # z stats
        max_z = max(z_list) if z_list else None
        top5_mean_z = topk_mean(z_list, 5) if z_list else None

        # specificity / concentration
        if deg > 0:
            max_species_fraction = max(neigh_species.values()) / deg
        else:
            max_species_fraction = 0.0

        H = entropy_from_counts(neigh_species)  # natural log
        # normalize entropy by log(n_species) (only defined if n_species >= 2)
        if n_species >= 2:
            H_norm = H / math.log(n_species)
        else:
            H_norm = 0.0

        n_eff = math.exp(H)  # effective species count
        part = participation_coefficient(neigh_species, deg)

        # clustering + betweenness
        cc = clustering.get(u) if clustering is not None else None
        bw = betweenness.get(u) if betweenness is not None else None

        feats.append(
            NodeFeatures(
                u=u,
                species=sp_u,
                component_id=component_id_of.get(u, -1),

                deg_xsp=deg,
                n_species=n_species,
                top1_jac=top1_jac,
                top5_mean_jac=top5_mean_jac,
                top1_shared=top1_shared,
                top5_mean_shared=top5_mean_shared,

                max_z=max_z,
                top5_mean_z=top5_mean_z,

                max_species_fraction=max_species_fraction,
                species_entropy=H,
                species_entropy_norm=H_norm,
                n_eff_species=n_eff,
                participation_coeff=part,

                clustering_coeff=cc,
                betweenness=bw,
            )
        )

    return feats


def write_node_features(path: Path, feats: List[NodeFeatures]) -> None:
    """Write protein_features.tsv."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "u", "species", "component_id",
            "deg_xsp", "n_species",
            "top1_jac", "top5_mean_jac",
            "top1_shared", "top5_mean_shared",
            "max_z", "top5_mean_z",
            "max_species_fraction",
            "species_entropy", "species_entropy_norm",
            "n_eff_species", "participation_coeff",
            "clustering_coeff", "betweenness",
        ])
        for x in feats:
            w.writerow([
                x.u, x.species, x.component_id,
                x.deg_xsp, x.n_species,
                f"{x.top1_jac:.6f}", f"{x.top5_mean_jac:.6f}",
                x.top1_shared, f"{x.top5_mean_shared:.6f}",
                "" if x.max_z is None else f"{x.max_z:.6f}",
                "" if x.top5_mean_z is None else f"{x.top5_mean_z:.6f}",
                f"{x.max_species_fraction:.6f}",
                f"{x.species_entropy:.6f}", f"{x.species_entropy_norm:.6f}",
                f"{x.n_eff_species:.6f}",
                f"{x.participation_coeff:.6f}",
                "" if x.clustering_coeff is None else f"{x.clustering_coeff:.6f}",
                "" if x.betweenness is None else f"{x.betweenness:.6f}",
            ])


# =========================
# Connected components + component-level HGT signals (addition #2)
# =========================

def compute_components(G) -> Tuple[Dict[str, int], List[List[str]]]:
    """
    Compute connected components of an undirected graph.

    Returns:
      component_id_of: dict node -> component_id (0..k-1)
      components: list of components, each a list of nodes

    Notes:
      - Components are returned in descending size order (largest first), for convenience.
      - component_id values follow that order.
    """
    if nx is None:
        raise RuntimeError("networkx is required for compute_components()")

    comps = [list(c) for c in nx.connected_components(G)]
    comps.sort(key=len, reverse=True)

    component_id_of: Dict[str, int] = {}
    for cid, nodes in enumerate(comps):
        for u in nodes:
            component_id_of[u] = cid

    return component_id_of, comps


def compute_component_features(
    G,
    components: List[List[str]],
    z_attr: str = "z_robust",
    z0: float = 3.0,
) -> List[ComponentFeatures]:
    """
    Compute component-level features:

    Species mixing:
      - n_species(C)
      - species_entropy(C) = -Σ p log p
      - species_entropy_norm(C) = H / log(n_species) (0 if n_species<2)

    High-z concentration (graph-level HGT signal):
      - high_z_frac(C) = fraction of edges in component with z >= z0
      - max_z(C)
      - top10_mean_z(C) over edges

    If a component has no numeric z values at all, the high-z fields are set to None.
    """
    feats: List[ComponentFeatures] = []

    for cid, nodes in enumerate(components):
        node_set = set(nodes)
        size = len(nodes)

        # ---- Species distribution and entropy ----
        species_counts = Counter(G.nodes[u].get("species", "") for u in nodes)
        n_species = sum(1 for s, c in species_counts.items() if c > 0)

        H = entropy_from_counts(species_counts)  # natural log entropy
        if n_species >= 2:
            H_norm = H / math.log(n_species)
        else:
            H_norm = 0.0

        # ---- Edge z statistics inside component ----
        num_edges = 0
        z_vals: List[float] = []
        high_z_count = 0

        # Iterate edges internal to the component efficiently:
        # Go over adjacency but count each undirected edge once (u < v).
        for u in nodes:
            for v, edata in G[u].items():
                if v not in node_set:
                    continue
                if u >= v:
                    continue  # ensure each undirected edge counted once
                num_edges += 1
                z = edata.get(z_attr)
                if z is None:
                    continue
                zf = float(z)
                z_vals.append(zf)
                if zf >= z0:
                    high_z_count += 1

        if z_vals:
            high_z_frac = high_z_count / num_edges if num_edges > 0 else 0.0
            max_z = max(z_vals)
            top10_mean_z = topk_mean(z_vals, 10)
        else:
            high_z_frac = None
            max_z = None
            top10_mean_z = None

        feats.append(
            ComponentFeatures(
                component_id=cid,
                size=size,
                num_edges=num_edges,
                n_species=n_species,
                species_entropy=H,
                species_entropy_norm=H_norm,
                high_z_frac=high_z_frac,
                max_z=max_z,
                top10_mean_z=top10_mean_z,
            )
        )

    return feats


def write_component_features(path: Path, feats: List[ComponentFeatures]) -> None:
    """Write component_features.tsv."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "component_id", "size", "num_edges",
            "n_species", "species_entropy", "species_entropy_norm",
            "high_z_frac", "max_z", "top10_mean_z",
        ])
        for c in feats:
            w.writerow([
                c.component_id, c.size, c.num_edges,
                c.n_species, f"{c.species_entropy:.6f}", f"{c.species_entropy_norm:.6f}",
                "" if c.high_z_frac is None else f"{c.high_z_frac:.6f}",
                "" if c.max_z is None else f"{c.max_z:.6f}",
                "" if c.top10_mean_z is None else f"{c.top10_mean_z:.6f}",
            ])


# =========================
# Candidate ranking / outputs
# =========================
def _minmax_scale(values: List[float]) -> List[float]:
    """Min-max scale to [0,1]. If constant, returns all zeros."""
    if not values:
        return []
    vmin = min(values)
    vmax = max(values)
    if vmax <= vmin:
        return [0.0] * len(values)
    return [(v - vmin) / (vmax - vmin) for v in values]

def score_hgt_likeness(
    node_feats: List[NodeFeatures],
    comp_feats: List[ComponentFeatures],
    min_component_size: int = 10,
    use_betweenness: bool = True,

    # --- Step-3 knobs ---
    H_min: float = 0.75,              # component entropy gate
    f_min: Optional[float] = None,    # optional node specificity gate (e.g. 0.4); None disables
    density_eta: float = 2.0,         # penalize dense components: (1-density)^eta
    cluster_eta: float = 1.0,         # reward bridge-like nodes: (1-clustering)^eta

    gamma_highz: float = 1.0,         # weight high_z_frac in context
    lambda_betweenness: float = 1.0,
    z_cap: float = 10.0,
    scale_floor_frac: float = 1e-3,
    context_floor: float = 1e-6,

    # --- Safety gates ---
    min_deg: int = 4,                 # reject leaves / tiny-degree nodes
    min_n_species: int = 2,           # require at least 2 neighbor species (local mixing)
    require_positive_z: bool = True,  # require top5_mean_z > 0 to have evidence
    entropy_node_floor: float = 0.15, # downweight near-zero neighbor-entropy (species-pure)
    require_any_anomaly: bool = True, # require Zz or Zf or Zb positive

    # --- NEW: evidence-strength + bridge-or-multispecies gates ---
    z_min: float = 2.0,               # require top5_mean_z >= z_min (nontrivial surprise)
    require_multispecies_or_bridge: bool = True,
    min_species_for_multispecies: int = 3,
    bw_min: float = 1e-12,            # treat bw > bw_min as "nonzero"
) -> List[Tuple[str, float]]:
    """
    Step-3 score (log-space, stable) with additional gates to prevent
    high-scoring weak artifacts (low z, weak similarity, leaf-like nodes).

    New gates:
      - require top5_mean_z >= z_min (default 1.0)
      - require (n_species >= 3) OR (betweenness > bw_min)  [if enabled]

    This should suppress:
      - nodes with tiny z like 0.02
      - nodes with only 2-species neighborhood and no bridging role
    """
    comp_by_id: Dict[int, ComponentFeatures] = {c.component_id: c for c in comp_feats}

    by_c: Dict[int, List[NodeFeatures]] = defaultdict(list)
    for nf in node_feats:
        by_c[nf.component_id].append(nf)

    def robust_center_scale(xs: List[float]) -> Tuple[float, float]:
        m = median(xs)
        d = mad(xs, m)
        s = 1.4826 * d
        floor = scale_floor_frac * (abs(m) + 1.0)
        return m, max(s, floor)

    comp_baselines: Dict[int, Dict[str, Tuple[float, float]]] = {}
    for cid, lst in by_c.items():
        z_vals = [(nf.top5_mean_z or 0.0) for nf in lst]
        f_vals = [nf.max_species_fraction for nf in lst]
        d_vals = [float(nf.deg_xsp) for nf in lst]
        if use_betweenness:
            b_vals = [(nf.betweenness or 0.0) for nf in lst]
        else:
            b_vals = [0.0 for _ in lst]

        comp_baselines[cid] = {
            "z": robust_center_scale(z_vals),
            "f": robust_center_scale(f_vals),
            "d": robust_center_scale(d_vals),
            "b": robust_center_scale(b_vals),
        }

    def clamp_pos(z: float) -> float:
        if z <= 0.0:
            return 0.0
        return z_cap if z > z_cap else z

    ranked: List[Tuple[str, float]] = []

    # log-space weights
    w_z = 1.0
    w_f = 2.0
    w_d = 0.5
    w_b = 1.0

    for nf in node_feats:
        cid = nf.component_id
        comp = comp_by_id.get(cid)
        if comp is None:
            ranked.append((nf.u, 0.0))
            continue

        comp_size = comp.size
        if comp_size < min_component_size:
            ranked.append((nf.u, 0.0))
            continue

        # Component entropy gate
        Hc = comp.species_entropy_norm
        if Hc < H_min:
            ranked.append((nf.u, 0.0))
            continue

        # Reject tiny-degree / low local mixing
        if nf.deg_xsp < min_deg:
            ranked.append((nf.u, 0.0))
            continue

        if nf.n_species < min_n_species:
            ranked.append((nf.u, 0.0))
            continue

        # Optional node specificity gate
        if f_min is not None and nf.max_species_fraction < f_min:
            ranked.append((nf.u, 0.0))
            continue

        # Require positive z evidence
        z_u_raw = nf.top5_mean_z
        if require_positive_z:
            if z_u_raw is None or z_u_raw <= 0.0:
                ranked.append((nf.u, 0.0))
                continue

        # NEW: require nontrivial surprise strength
        if z_u_raw is None or z_u_raw < z_min:
            ranked.append((nf.u, 0.0))
            continue

        # NEW: require multi-species neighborhood OR bridge evidence
        if require_multispecies_or_bridge:
            bw = (nf.betweenness or 0.0) if use_betweenness else 0.0
            if nf.n_species < min_species_for_multispecies and bw <= bw_min:
                ranked.append((nf.u, 0.0))
                continue

        # Component density penalty
        if comp_size >= 2:
            dens = (2.0 * comp.num_edges) / (comp_size * (comp_size - 1.0))
        else:
            dens = 0.0
        dens = min(1.0, max(0.0, dens))
        dens_pen = max(1e-6, (1.0 - dens) ** density_eta)

        # Context: emphasize high-z concentration
        highz = comp.high_z_frac if comp.high_z_frac is not None else 0.0
        context = max(highz ** gamma_highz, context_floor)

        # Evidence gates
        size_gate = math.log1p(comp_size)
        deg_gate = math.log1p(nf.deg_xsp)
        if deg_gate <= 0.0:
            ranked.append((nf.u, 0.0))
            continue

        # Bridge proxy from clustering (low clustering => bridge-ish)
        cc = nf.clustering_coeff if nf.clustering_coeff is not None else 1.0
        cc = min(1.0, max(0.0, cc))
        bridge_local = max(1e-6, (1.0 - cc) ** cluster_eta)

        # Downweight very low neighbor-entropy (species-pure)
        ent_u = nf.species_entropy_norm
        if ent_u < entropy_node_floor:
            ent_pen = 0.1 + 0.9 * (ent_u / max(entropy_node_floor, 1e-9))
        else:
            ent_pen = 1.0

        # Component-relative deviations
        z_u = nf.top5_mean_z or 0.0
        f_u = nf.max_species_fraction
        d_u = float(nf.deg_xsp)
        b_u = (nf.betweenness or 0.0) if use_betweenness else 0.0

        (mz, sz) = comp_baselines[cid]["z"]
        (mf, sf) = comp_baselines[cid]["f"]
        (md, sd) = comp_baselines[cid]["d"]
        (mb, sb) = comp_baselines[cid]["b"]

        Zz = clamp_pos((z_u - mz) / sz)
        Zf = clamp_pos((f_u - mf) / sf)
        Zd = clamp_pos((d_u - md) / sd)
        Zb = clamp_pos((b_u - mb) / sb) if use_betweenness else 0.0

        # Require at least one anomaly relative to component
        if require_any_anomaly and (Zz == 0.0 and Zf == 0.0 and Zb == 0.0):
            ranked.append((nf.u, 0.0))
            continue

        # log-space aggregation
        logS = 0.0
        logS += math.log(size_gate)
        logS += math.log(context)
        logS += math.log(deg_gate)
        logS += math.log(dens_pen)
        logS += math.log(bridge_local)
        logS += math.log(ent_pen)

        logS += w_z * math.log1p(Zz)
        logS += w_f * math.log1p(Zf)
        logS += w_d * math.log1p(0.25 * Zd)
        if use_betweenness:
            logS += w_b * math.log1p(lambda_betweenness * Zb)

        ranked.append((nf.u, math.exp(logS)))

    ranked.sort(key=lambda t: t[1], reverse=True)
    return ranked




def write_hgt_candidates(
    path: Path,
    node_feats: List[NodeFeatures],
    scores: Dict[str, float],
    top_n: int = 200,
) -> None:
    """
    Write hgt_candidates.tsv containing top N proteins with key features.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # Sort nodes by score
    rows = sorted(
        node_feats,
        key=lambda nf: scores.get(nf.u, 0.0),
        reverse=True
    )[:top_n]

    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "u", "species", "component_id", "score",
            "deg_xsp", "n_species",
            "max_z", "top5_mean_z",
            "max_species_fraction", "species_entropy_norm",
            "betweenness", "clustering_coeff",
        ])

        for nf in rows:
            score = scores.get(nf.u, 0.0)

            w.writerow([
                nf.u,
                nf.species,
                nf.component_id,
                f"{score:.6f}",
                nf.deg_xsp,
                nf.n_species,
                "" if nf.max_z is None else f"{nf.max_z:.6f}",
                "" if nf.top5_mean_z is None else f"{nf.top5_mean_z:.6f}",
                f"{nf.max_species_fraction:.6f}",
                f"{nf.species_entropy_norm:.6f}",
                "" if nf.betweenness is None else f"{nf.betweenness:.6f}",
                "" if nf.clustering_coeff is None else f"{nf.clustering_coeff:.6f}",
            ])


# =========================
# Orchestration (main pipeline)
# =========================
def run_pipeline(
    in_path: Path,
    out_dir: Path,
    weight_for_z: str = "jaccard",
    min_pair_edges_for_z: int = 30,
    z0: float = 3.0,
    compute_betweenness: bool = True,
) -> None:
    """
    Run the full analysis pipeline with per-stage logging + timing.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    with stage("1) Read edges"):
        edges = read_edges(in_path)
        edges = dedupe_edges(edges)
        log_counts(edges=edges)

    with stage("2) Build graph"):
        G = build_graph(edges, use_networkx=True)
        log_counts(G=G)

    with stage("3) Robust per-species-pair stats"):
        pair_stats = compute_pair_robust_stats(
            edges,
            weight=weight_for_z,
            min_pair_edges=min_pair_edges_for_z,
        )
        log_counts(pair_stats=pair_stats)
        print_pair_stats_sanity_table(
            pair_stats,
            min_pair_edges_for_z=min_pair_edges_for_z,
            top_n=15,
        )

    with stage("4) Edge features (z) + write edge_features.tsv"):
        edge_feats = compute_edge_features(
            edges,
            pair_stats=pair_stats,
            weight=weight_for_z,
            min_pair_edges_for_z=min_pair_edges_for_z,
        )
        write_edge_features(out_dir / "edge_features.tsv", edge_feats)
        print(f"[INFO] wrote {out_dir / 'edge_features.tsv'}")

    with stage("5) Attach z to graph edges"):
        attach_edge_z_to_graph(G, edge_feats)

    # sanity: count how many edges have z attached
    n_total = G.number_of_edges()
    n_with_z = sum(
        1 for _, _, d in G.edges(data=True)
        if "z_robust" in d
    )
    print(f"[SANITY] edges with z: {n_with_z}/{n_total}")

    with stage("6) Connected components"):
        component_id_of, components = compute_components(G)
        print(f"[INFO] components={len(components)}")

    with stage("7) Component features + write component_features.tsv"):
        comp_feats = compute_component_features(G, components, z0=z0)
        write_component_features(out_dir / "component_features.tsv", comp_feats)
        print(f"[INFO] wrote {out_dir / 'component_features.tsv'}")

    sizes = [len(c) for c in components]
    print(
        f"[SANITY] #components={len(components)} | largest={sizes[0]} | median={sizes[len(sizes) // 2] if sizes else 0}")
    mix = [c for c in comp_feats if c.n_species >= 3 and c.size >= 20]
    mix.sort(key=lambda x: (x.species_entropy_norm, x.size), reverse=True)
    print(f"[SANITY] mixed comps (>=3 species, >=20 nodes): {len(mix)}")
    for c in mix[:5]:
        print(
            f"  cid={c.component_id} size={c.size} n_species={c.n_species} Hn={c.species_entropy_norm:.3f} highz={'' if c.high_z_frac is None else f'{c.high_z_frac:.3f}'}")

    with stage("8) Node features + write protein_features.tsv"):
        node_feats = compute_node_features(
            G,
            component_id_of=component_id_of,
            compute_betweenness=compute_betweenness,
            betweenness_mode="high_z",
            z_threshold_for_high_subgraph=z0,
        )
        write_node_features(out_dir / "protein_features.tsv", node_feats)
        print(f"[INFO] wrote {out_dir / 'protein_features.tsv'}")

    with stage("9) Rank proteins + write hgt_candidates.tsv"):
        ranked = score_hgt_likeness(node_feats, comp_feats, min_component_size=10,
                                                                use_betweenness=compute_betweenness)
        scores = {u: s for (u, s) in ranked}
        write_hgt_candidates(out_dir / "hgt_candidates.tsv", node_feats, scores, top_n=200)
        print(f"[INFO] wrote {out_dir / 'hgt_candidates.tsv'}")

    print("[SANITY] top 10 HGT-like proteins:")
    for u, s in ranked[:10]:
        print(f"  {s:10.3f}  {u}")

    print("[SANITY] score stats:",
          "max=", ranked[0][1],
          "p50=", ranked[len(ranked) // 2][1],
          "p95=", ranked[int(0.95 * len(ranked))][1])

    vals = [s for _, s in ranked if s > 0]
    vals.sort()
    print(
        f"[SANITY] score stats: n={len(vals)} min={vals[0]:.6g} p50={vals[len(vals) // 2]:.6g} p95={vals[int(0.95 * len(vals))]:.6g} max={vals[-1]:.6g}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Graph construction + HGT-like analysis (stubs).")
    ap.add_argument("--in_edges", required=True, type=Path, help="Input CSV/TSV of undirected pruned edges")
    ap.add_argument("--out_dir", required=True, type=Path, help="Output directory for TSV deliverables")
    ap.add_argument("--weight_for_z", choices=["jaccard", "shared"], default="jaccard")
    ap.add_argument("--min_pair_edges_for_z", type=int, default=30)
    ap.add_argument("--z0", type=float, default=3.0, help="High-z threshold for component concentration + subgraph")
    ap.add_argument("--no_betweenness", action="store_true", help="Skip betweenness (fast)")
    args = ap.parse_args()

    run_pipeline(
        in_path=args.in_edges,
        out_dir=args.out_dir,
        weight_for_z=args.weight_for_z,
        min_pair_edges_for_z=args.min_pair_edges_for_z,
        z0=args.z0,
        compute_betweenness=not args.no_betweenness,
    )


if __name__ == "__main__":
    main()
