import matplotlib.pyplot as plt
from pathlib import Path
from Bio import Phylo
from typing import Dict, List, Tuple
from matplotlib.collections import LineCollection
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from ..constants import PHYLO_REG_FONT_SIZE, PHYLO_TITLE


def create_dist_mat(n: int, names: List[str], pairwise_identity: Dict[Tuple[str, str], float]) -> List[List[float]]:
    """
    Create a distance matrix from pairwise identity data.
    :param n: Number of names.
    :param names: List of names.
    :param pairwise_identity: Dictionary mapping (name1, name2) to their pairwise identity.
    :return: A lower-triangular distance matrix.
    """
    mat: List[List[float]] = []
    for i in range(n):
        row: List[float] = []
        for j in range(i + 1):
            if i == j:
                dist = 0.0
            else:
                a, b = names[i], names[j]
                ident = pairwise_identity.get((a, b), pairwise_identity.get((b, a), None))
                if ident is None:
                    dist = 1.0
                else:
                    dist = 1.0 - float(ident)
                    dist = max(0.0, min(1.0, dist))
            row.append(dist)
        mat.append(row)
    return mat


def export_nj_tree(node_ids: List[str], node_to_organism: Dict[str, str], pairwise_identity: Dict[Tuple[str, str], float],
                   out_png: str, title: str = PHYLO_TITLE) -> None:
    """
    Create and export a neighbor-joining tree based on pairwise identity.
    :param node_ids: The list of node IDs to include in the tree.
    :param node_to_organism: A mapping from node IDs to organism names.
    :param pairwise_identity: A dictionary mapping (node_id1, node_id2) to their pairwise identity (0.0 to 1.0).
    :param out_png: Path to save the output PNG image of the tree.
    :param title: Title for the tree plot.
    """
    names: List[str] = list(node_ids)
    n = len(names)
    if n < 2:
        return

    tri: List[List[float]] = create_dist_mat(n, names, pairwise_identity)

    dm = DistanceMatrix(names=names, matrix=tri)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    for clade in tree.get_terminals():
        prot_id = clade.name
        org = node_to_organism.get(prot_id, prot_id)
        clade.name = org

    plt.rcParams['path.simplify'] = False
    plt.rcParams['font.family'] = 'Calibri'

    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title(title, fontsize=24, fontweight='bold', pad=12)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='both', width=1.5, length=6)
    ax.set_xlabel(ax.get_xlabel(), fontsize=PHYLO_REG_FONT_SIZE, fontweight='bold')
    ax.set_ylabel(ax.get_ylabel(), fontsize=PHYLO_REG_FONT_SIZE, fontweight='bold')

    Phylo.draw(tree, do_show=False, axes=ax)

    for text in ax.texts:
        if not text.get_text().strip().startswith('Inner'):  # Skip internal node labels
            text.set_fontsize(PHYLO_REG_FONT_SIZE)
            text.set_fontweight('bold')

    for artist in ax.get_children():
        if isinstance(artist, LineCollection):
            artist.set_linewidth(4.0)
            artist.set_antialiased(True)

    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)