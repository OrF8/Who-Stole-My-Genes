from collections import Counter, defaultdict
from typing import Dict


def add_name_idx(node_to_organism: Dict[str, str]) -> None:
    """
    Add an index to organism names to ensure uniqueness.
    :param node_to_organism: A mapping from node IDs to organism names.
    """
    frequencies = Counter(node_to_organism.values())
    counters: Dict[str, int] = defaultdict(int)

    for node, org in node_to_organism.items():
        if frequencies[org] > 1:
            counters[org] += 1
            node_to_organism[node] = f"{org} {counters[org]}"
        else:
            node_to_organism[node] = org


def normalize_species_name(organism: str) -> str:
    """
    Normalize organism name to "Genus species" format.
    :param organism: The input organism name.
    :return: The normalized organism name.
    """
    parts = organism.strip().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else organism.strip()