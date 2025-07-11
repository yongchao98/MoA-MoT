import collections

def count_3_matchings(num_vertices, adjacency_list):
    """
    Counts the number of 3-matchings in a graph.

    This function implements a polynomial-time algorithm based on iterating
    through each edge and counting 2-matchings in the remaining graph.
    The complexity is O(num_edges * num_vertices), which is O(n^3) for dense graphs.

    A 3-matching is a set of 3 edges with no shared vertices.
    The method is based on the formula:
    3 * M_3 = sum_{(u,v) in E} M_2(G - {u,v})
    where M_k(G) is the number of k-matchings in graph G, and G - {u,v}
    is the graph with vertices u, v and all incident edges removed.

    The number of 2-matchings M_2 in a graph H is given by:
    M_2(H) = (|E(H)| choose 2) - sum_{w in V(H)} (degree_H(w) choose 2)

    Args:
        num_vertices (int): The number of vertices in the graph, V = {0, 1, ..., n-1}.
        adjacency_list (list of lists): The adjacency list representation of the graph.
                                        adjacency_list[i] contains neighbors of vertex i.

    Returns:
        None: This function prints the result instead of returning it.
    """
    if num_vertices < 6:
        print("Graph has fewer than 6 vertices, so no 3-matching is possible.")
        print(f"Final count of 3-matchings: 0")
        return

    # For faster neighbor lookups, convert lists to sets
    adj_sets = [set(neighbors) for neighbors in adjacency_list]
    degrees = [len(neighbors) for neighbors in adjacency_list]
    
    edges = []
    for i in range(num_vertices):
        for j in adjacency_list[i]:
            if i < j:
                edges.append((i, j))
    
    num_edges = len(edges)
    
    # This will accumulate the count of 3-matchings multiplied by 3
    total_3_matchings_times_3 = 0
    
    # Iterate over each edge e1 = (u, v)
    for u, v in edges:
        # Now, count 2-matchings in the graph G' = G - {u,v}.
        
        # Calculate the number of edges in G'.
        # Start with total edges and subtract edges connected to u or v.
        # An edge is counted twice (deg(u)+deg(v)), but (u,v) itself is counted
        # once for each endpoint, so we add 1 back since it was subtracted twice.
        num_edges_in_guv = num_edges - degrees[u] - degrees[v] + 1
        
        # This is the C(|E(G')|, 2) term
        term1 = (num_edges_in_guv * (num_edges_in_guv - 1)) // 2
        
        # This will be sum_{w} C(deg_G'(w), 2)
        term2 = 0
        
        # Iterate over all vertices w in G'
        for w in range(num_vertices):
            if w == u or w == v:
                continue
            
            # Calculate degree of w in G'
            deg_w_in_guv = degrees[w]
            if u in adj_sets[w]:
                deg_w_in_guv -= 1
            if v in adj_sets[w]:
                deg_w_in_guv -= 1
            
            term2 += (deg_w_in_guv * (deg_w_in_guv - 1)) // 2
        
        # Number of 2-matchings in G'
        num_2_matchings_in_guv = term1 - term2
        total_3_matchings_times_3 += num_2_matchings_in_guv
        
    # Each 3-matching {e1, e2, e3} is counted 3 times:
    # once when e1 is chosen, once for e2, and once for e3.
    # So we divide the total sum by 3.
    num_3_matchings = total_3_matchings_times_3 // 3
    
    print(f"The calculation for the number of 3-matchings (M_3):")
    print(f"Let S be the sum of 2-matchings over all edge-deleted subgraphs, S = sum_{{(u,v) in E}} M_2(G-{{u,v}})")
    print(f"We found S = {total_3_matchings_times_3}")
    print(f"The number of 3-matchings is given by the equation M_3 = S / 3.")
    print(f"{total_3_matchings_times_3} / 3 = {num_3_matchings}")
    print(f"\nFinal count of 3-matchings: {num_3_matchings}")


# --- Example Usage ---
# We will test the function on a simple graph.
# Consider a graph consisting of a 6-vertex cycle: 0-1-2-3-4-5-0.
# Edges: (0,1), (1,2), (2,3), (3,4), (4,5), (5,0)
# This graph has two 3-matchings: {(0,1), (2,3), (4,5)} and {(1,2), (3,4), (5,0)}.
n_vertices_example = 6
adj_list_example = [
    [1, 5],
    [0, 2],
    [1, 3],
    [2, 4],
    [3, 5],
    [0, 4]
]

print("--- Counting 3-matchings in a 6-cycle graph ---")
count_3_matchings(n_vertices_example, adj_list_example)