import networkx as nx
from itertools import combinations

def count_k_cliques(G, k):
    """Counts the number of k-cliques in G."""
    if k < 2:
        return len(G.nodes()) if k == 1 else 0
    
    count = 0
    nodes = list(G.nodes())
    for subset in combinations(nodes, k):
        subgraph = G.subgraph(subset)
        if subgraph.number_of_edges() == k * (k - 1) // 2:
            count += 1
    return count

def count_induced_k_matchings(G, k):
    """Counts the number of induced k-matchings in G."""
    if k == 0:
        return 1
    if k < 1:
        return 0

    count = 0
    edges = list(G.edges())
    
    # Iterate through all combinations of k edges
    for edge_subset in combinations(edges, k):
        vertex_set = set()
        is_matching = True
        
        # Check if it is a matching
        for u, v in edge_subset:
            if u in vertex_set or v in vertex_set:
                is_matching = False
                break
            vertex_set.add(u)
            vertex_set.add(v)
        
        if not is_matching:
            continue
            
        # Check if the matching is induced
        induced_subgraph = G.subgraph(vertex_set)
        if induced_subgraph.number_of_edges() == k:
            count += 1
            
    return count

def count_induced_k_k_bicliques(G, k):
    """Counts the number of induced k-by-k-bicliques in G."""
    if k < 1:
        return 0

    count = 0
    nodes = list(G.nodes())
    
    # Iterate through all disjoint pairs of k-sized vertex sets
    for a_nodes in combinations(nodes, k):
        remaining_nodes = [n for n in nodes if n not in a_nodes]
        if len(remaining_nodes) < k:
            continue
            
        for b_nodes in combinations(remaining_nodes, k):
            a_set = set(a_nodes)
            b_set = set(b_nodes)
            vertex_set = a_set.union(b_set)
            
            subgraph = G.subgraph(vertex_set)
            
            # Check induced biclique conditions
            is_biclique = True
            
            # 1. No edges within A
            if subgraph.subgraph(a_set).number_of_edges() > 0:
                is_biclique = False
            
            # 2. No edges within B
            if is_biclique and subgraph.subgraph(b_set).number_of_edges() > 0:
                is_biclique = False

            # 3. All edges between A and B must exist
            if is_biclique:
                expected_edges = k * k
                actual_edges = 0
                for u in a_set:
                    for v in b_set:
                        if G.has_edge(u, v):
                            actual_edges += 1
                if actual_edges != expected_edges:
                    is_biclique = False
            
            if is_biclique:
                count += 1
    return count

def solve_pcount(G, k):
    """
    Solves the PCount problem for a given graph G and parameter k.
    This demonstrates a non-FPT (brute-force) approach.
    """
    print(f"Solving PCount for k={k} on a graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    num_cliques = count_k_cliques(G, k)
    num_induced_matchings = count_induced_k_matchings(G, k)
    num_induced_bicliques = count_induced_k_k_bicliques(G, k)
    
    total = num_cliques + num_induced_matchings + num_induced_bicliques
    
    # Output the final equation with each number
    print("Final Equation:")
    print(f"{num_cliques} (k-cliques) + {num_induced_matchings} (induced k-matchings) + {num_induced_bicliques} (induced k-by-k-bicliques) = {total}")
    
    return total

if __name__ == '__main__':
    # Example Usage: A cycle graph on 5 vertices (C5) and k=2
    # This graph has five 2-cliques (edges), one induced 2-matching, and zero induced 2x2-bicliques.
    # Expected: 5 + 1 + 0 = 6
    G_example = nx.cycle_graph(5)
    k_example = 2
    
    solve_pcount(G_example, k_example)

    print("\n" + "="*40 + "\n")
    
    # Another example: A complete graph on 4 vertices (K4) and k=3
    # This graph has four 3-cliques, zero induced 3-matchings (needs 6 vertices), and zero induced 3x3-bicliques.
    # Expected: 4 + 0 + 0 = 4
    G_example_2 = nx.complete_graph(4)
    k_example_2 = 3
    solve_pcount(G_example_2, k_example_2)
