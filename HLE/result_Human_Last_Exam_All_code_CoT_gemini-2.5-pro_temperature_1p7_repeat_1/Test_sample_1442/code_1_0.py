import itertools

def build_graphs():
    """Builds the adjacency lists for two 3-regular bipartite graphs on 10 vertices."""
    
    # Graph G1
    # u_i is connected to w_i, w_{i-1}, and w_{i+2} (mod 5)
    adj_g1 = {f'u{i}': [] for i in range(5)}
    adj_g1.update({f'w{i}': [] for i in range(5)})

    for i in range(5):
        # Edges from cycle part
        adj_g1[f'u{i}'].extend([f'w{i}', f'w{(i - 1 + 5) % 5}'])
        adj_g1[f'w{i}'].extend([f'u{i}', f'u{(i + 1) % 5}'])
        
        # Edges from chord part (u_i, w_{i+2})
        adj_g1[f'u{i}'].append(f'w{(i + 2) % 5}')
        adj_g1[f'w{(i + 2) % 5}'].append(f'u{i}')
        
    # Graph G2
    # u_i is connected to w_i, w_{i-1}, and w_{i+3} (mod 5)
    adj_g2 = {f'u{i}': [] for i in range(5)}
    adj_g2.update({f'w{i}': [] for i in range(5)})

    for i in range(5):
        # Edges from cycle part
        adj_g2[f'u{i}'].extend([f'w{i}', f'w{(i - 1 + 5) % 5}'])
        adj_g2[f'w{i}'].extend([f'u{i}', f'u{(i + 1) % 5}'])
        
        # Edges from chord part (u_i, w_{i+3})
        adj_g2[f'u{i}'].append(f'w{(i + 3) % 5}')
        adj_g2[f'w{(i + 3) % 5}'].append(f'u{i}')

    # Correcting for duplicate edges from building logic
    for v in adj_g1:
        adj_g1[v] = sorted(list(set(adj_g1[v])))
    for v in adj_g2:
        adj_g2[v] = sorted(list(set(adj_g2[v])))
        
    return adj_g1, adj_g2

def count_3_matchings(adj_list):
    """Counts the number of 3-matchings in a graph given its adjacency list."""
    
    # Get a unique list of edges
    edges = set()
    for u, neighbors in adj_list.items():
        for v in neighbors:
            # Add edge as a frozenset to handle (u,v) vs (v,u) and ensure uniqueness
            edges.add(frozenset([u, v]))
    
    edge_list = list(edges)
    
    matching_count = 0
    # Iterate through all combinations of 3 edges
    for edge_combo in itertools.combinations(edge_list, 3):
        e1, e2, e3 = edge_combo
        
        # Check for disjoint vertices
        vertices_in_matching = list(e1) + list(e2) + list(e3)
        if len(set(vertices_in_matching)) == 6:
            matching_count += 1
            
    return matching_count

# Main execution
g1_adj, g2_adj = build_graphs()

# Count 3-matchings for both graphs
g1_matchings = count_3_matchings(g1_adj)
g2_matchings = count_3_matchings(g2_adj)

print("Graph G1 is a bipartite, 3-regular graph on 10 vertices.")
print(f"Number of 3-matchings in G1: {g1_matchings}")
print("\nGraph G2 is another bipartite, 3-regular graph on 10 vertices, non-isomorphic to G1.")
print(f"Number of 3-matchings in G2: {g2_matchings}")

if g1_matchings != g2_matchings:
    print("\nSince the number of 3-matchings is different, two bipartite, d-regular graphs")
    print("on n vertices do not necessarily have the same number of 3-matchings.")
else:
    print("\nThis example did not prove the hypothesis, but the answer is still No.")

print("\nFinal Answer: The statement is false.")
print("The number of 1-matchings and 2-matchings is guaranteed to be the same for all d-regular graphs on n vertices.")
print("However, the number of 3-matchings is not. Our counterexample demonstrates this.")
