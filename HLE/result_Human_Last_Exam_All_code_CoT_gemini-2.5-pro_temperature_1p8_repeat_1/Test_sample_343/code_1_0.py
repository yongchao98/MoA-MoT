import networkx as nx
from itertools import combinations

def count_cliques(G, k):
    """Counts the number of k-cliques in a graph G."""
    if k == 0:
        return 1
    if k > len(G.nodes):
        return 0
    
    count = 0
    for node_combo in combinations(G.nodes(), k):
        subgraph = G.subgraph(node_combo)
        if subgraph.number_of_edges() == k * (k - 1) // 2:
            count += 1
    return count

def count_induced_matchings(G, k):
    """Counts the number of induced k-matchings in a graph G."""
    if k == 0:
        return 1
    if 2 * k > len(G.nodes):
        return 0
        
    count = 0
    for node_combo in combinations(G.nodes(), 2 * k):
        subgraph = G.subgraph(node_combo)
        # An induced k-matching is a graph on 2k vertices where every vertex has degree 1
        # (i.e., a k-regular graph which is also a perfect matching)
        if all(d == 1 for _, d in subgraph.degree()):
            count += 1
    return count
    
def count_induced_bicliques(G, k):
    """Counts the number of induced k-by-k-bicliques in a graph G."""
    if k == 0:
        return 1
    if 2 * k > len(G.nodes):
        return 0

    count = 0
    for node_combo in combinations(G.nodes(), 2 * k):
        subgraph = G.subgraph(node_combo)
        # Check all partitions into two sets of size k
        nodes = list(node_combo)
        for part_a_indices in combinations(range(2 * k), k):
            part_a = {nodes[i] for i in part_a_indices}
            # We only need to check one representative partition (A, B) not (B, A)
            # The combination selection for A takes care of this if we start from a fixed node.
            # To avoid double counting A,B and B,A, we fix the first element of node_combo to be in part A.
            if 0 not in part_a_indices:
                continue

            part_b = set(nodes) - part_a

            is_biclique = True
            # Check for no edges within part A
            if G.subgraph(part_a).number_of_edges() > 0:
                is_biclique = False
            # Check for no edges within part B
            if is_biclique and G.subgraph(part_b).number_of_edges() > 0:
                is_biclique = False
            
            # Check for all edges between part A and part B
            if is_biclique:
                for u in part_a:
                    for v in part_b:
                        if not G.has_edge(u, v):
                            is_biclique = False
                            break
                    if not is_biclique:
                        break
            
            if is_biclique:
                count += 1
    return count

def p_count(G, k):
    """Calculates the PCount of a graph G for parameter k."""
    num_cliques = count_cliques(G, k)
    num_induced_matchings = count_induced_matchings(G, k)
    num_induced_bicliques = count_induced_bicliques(G, k)
    
    total = num_cliques + num_induced_matchings + num_induced_bicliques
    
    print(f"For k = {k}:")
    print(f"Number of {k}-cliques: {num_cliques}")
    print(f"Number of induced {k}-matchings: {num_induced_matchings}")
    print(f"Number of induced {k}-by-{k}-bicliques: {num_induced_bicliques}")
    print(f"Total PCount: {num_cliques} + {num_induced_matchings} + {num_induced_bicliques} = {total}")
    return total

if __name__ == '__main__':
    # Example Graph: A 6-cycle with a chord creating two triangles (K3)
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 5), (5, 3)])
    
    k = 3
    p_count(G, k)
