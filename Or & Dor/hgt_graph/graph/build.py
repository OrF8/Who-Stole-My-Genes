import networkx as nx
from typing import Dict, List, Tuple, Union

def build_nx_graph(node_attrs: Dict[str, Dict[str, Union[str, int]]],
                   edges: List[Tuple[str, str, Dict[str, float]]]) -> nx.Graph:
    """
    Build a NetworkX graph from node attributes and edges.
    :param node_attrs: A dictionary of node attributes.
    :param edges: A list of edges with weights.
    :return: A NetworkX graph.
    """
    G = nx.Graph()
    for node_id, attrs in node_attrs.items():
        G.add_node(node_id, **attrs)
    for u, v, eattrs in edges:
        G.add_edge(u, v, **eattrs)
    return G
