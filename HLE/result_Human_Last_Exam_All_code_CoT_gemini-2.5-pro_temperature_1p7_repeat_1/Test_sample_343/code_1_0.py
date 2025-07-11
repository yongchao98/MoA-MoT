import networkx as nx
from itertools import combinations

def count_k_cliques(G, k):
    """Counts the number of k-cliques in G."""
    count = 0
    # A k-clique is a set of k vertices where every two distinct vertices are adjacent.
    for nodes_subset in combinations(G.nodes(), k):
        is_clique = True
        # Check all pairs of vertices in the subset for an edge.
        for u, v in combinations(nodes_subset, 2):
            if not G.has_edge(u, v):
                is_clique = False
                break
        if is_clique:
            count += 1
    return count

def count_induced_k_matchings(G, k):
    """Counts the number of induced k-matchings in G."""
    count = 0
    edges = list(G.edges())
    if k == 0:
        return 1
    # Iterate over all subsets of k edges.
    for edge_subset in combinations(edges, k):
        vertex_set = set()
        is_matching = True
        # Check if the edges form a matching (no shared vertices).
        for u, v in edge_subset:
            if u in vertex_set or v in vertex_set:
                is_matching = False
                break
            vertex_set.add(u)
            vertex_set.add(v)
        
        if not is_matching:
            continue
            
        # Check if the matching is induced.
        # This means no extra edges exist between the 2k vertices.
        is_induced = True
        # Check all pairs of vertices from the matching's endpoints.
        for u, v in combinations(vertex_set, 2):
            # If an edge exists in G, it must be one of the chosen matching edges.
            if G.has_edge(u, v):
                if (u, v) not in edge_subset and (v, u) not in edge_subset:
                    is_induced = False
                    break
        if not is_induced:
            continue
            
        if is_induced:
            count += 1
            
    return count

def count_induced_k_by_k_bicliques(G, k):
    """Counts the number of induced k-by-k-bicliques (K_k,k) in G."""
    count = 0
    nodes = list(G.nodes())
    if k == 0:
        return 1
    if len(nodes) < 2 * k:
        return 0
    
    # Iterate over all ways to choose the first partition 'A' of size k.
    for nodes_A in combinations(nodes, k):
        A = set(nodes_A)
        remaining_nodes = [n for n in nodes if n not in A]
        if len(remaining_nodes) < k:
            continue
        
        # Iterate over all ways to choose the second partition 'B' of size k.
        for nodes_B in combinations(remaining_nodes, k):
            B = set(nodes_B)
            
            is_induced_biclique = True
            
            # 1. Check for no edges within A (A is an independent set).
            for u, v in combinations(A, 2):
                if G.has_edge(u, v):
                    is_induced_biclique = False
                    break
            if not is_induced_biclique: continue
            
            # 2. Check for no edges within B (B is an independent set).
            for u, v in combinations(B, 2):
                if G.has_edge(u, v):
                    is_induced_biclique = False
                    break
            if not is_induced_biclique: continue
            
            # 3. Check for all edges between A and B (A and B form a complete bipartite graph).
            for u in A:
                for v in B:
                    if not G.has_edge(u, v):
                        is_induced_biclique = False
                        break
                if not is_induced_biclique: break
            
            if is_induced_biclique:
                count += 1
    # Each K_k,k is found twice by this method (once as (A,B), once as (B,A)).
    # So we divide by 2.
    return count // 2 if k > 0 else 1

if __name__ == '__main__':
    # Let's create an example graph: a cycle of 5 vertices with one chord.
    # This graph is C5 + edge(0,2).
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2)])
    k = 3

    # Calculate the number of each substructure.
    num_cliques = count_k_cliques(G, k)
    num_induced_matchings = count_induced_k_matchings(G, k)
    num_induced_bicliques = count_induced_k_by_k_bicliques(G, k)
    
    # The total count for PCount.
    total_count = num_cliques + num_induced_matchings + num_induced_bicliques

    print(f"For k = {k}:")
    print(f"Number of k-cliques: {num_cliques}")
    print(f"Number of induced k-matchings: {num_induced_matchings}")
    print(f"Number of induced k-by-k-bicliques: {num_induced_bicliques}")
    print("---")
    # Output the final equation as requested
    print(f"{num_cliques} + {num_induced_matchings} + {num_induced_bicliques} = {total_count}")
