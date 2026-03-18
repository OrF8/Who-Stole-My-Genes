import os
import networkx as nx
from hgt_graph.taxonomy.cache import get_or_build_taxonomy_cache
from hgt_graph.taxonomy.annotate import taxonomic_distance, annotate_graph_with_taxonomy
from hgt_graph.constants import ORGANISM_KEY, NCBI_SLEEP
from typing import Dict, Tuple

if __name__ == '__main__':
    org_names_path: str = os.path.join('Helping Noam', os.path.join('data', 'organism_names.txt'))

    with open(org_names_path, 'r') as f:
        orgs = f.read().splitlines()

    node_attrs = {str(i): {ORGANISM_KEY: org } for i, org in enumerate(orgs)}

    taxonomy_cache = get_or_build_taxonomy_cache(
        cache_path=os.path.join('Helping Noam', os.path.join('data', 'taxonomy_cache.json')),
        node_attrs=node_attrs,
        save_path=os.path.join('Helping Noam', os.path.join('data', 'taxonomy_cache.json')),
        sleep_seconds=NCBI_SLEEP
    )

    G = nx.Graph()
    G.add_nodes_from(node_attrs.items())

    annotate_graph_with_taxonomy(
        G=G,
        taxonomy_cache=taxonomy_cache
    )

    distances: Dict[Tuple[str, str], int] = {}
    for i in range(len(orgs)):
        for j in range(i + 1, len(orgs)):
            u_id = str(i)
            v_id = str(j)
            d = taxonomic_distance(
                G.nodes[u_id],
                G.nodes[v_id]
            )
            distances[(orgs[i], orgs[j])] = d
    save_path: str = os.path.join('Helping Noam', 'tax_distances.txt')
    with open(save_path, 'w') as f:
        for (org1, org2), dist in distances.items():
            f.write(f'{org1},{org2},{dist}\n')

