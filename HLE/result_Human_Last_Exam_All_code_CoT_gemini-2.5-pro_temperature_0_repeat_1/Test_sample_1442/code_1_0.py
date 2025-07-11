import itertools

def count_3_matchings(graph_name, adj_list):
    """
    Calculates the number of 3-matchings in a graph.
    A 3-matching is a set of 3 edges with no shared vertices.
    """
    # Step 1: Create a list of unique edges from the adjacency list.
    # An edge is represented as a sorted tuple to avoid duplicates, e.g., (u, v) and (v, u).
    edges = set()
    for u, neighbors in enumerate(adj_list):
        for v in neighbors:
            edge = tuple(sorted((u, v)))
            edges.add(edge)
    
    edge_list = list(edges)
    
    # Step 2: Iterate through all combinations of 3 edges.
    k = 3
    matching_count = 0
    for edge_combo in itertools.combinations(edge_list, k):
        # Step 3: Check if the 3 edges are disjoint.
        # This is true if the 6 vertices involved are all unique.
        vertices_in_combo = set()
        for edge in edge_combo:
            vertices_in_combo.update(edge)
        
        if len(vertices_in_combo) == 2 * k:
            matching_count += 1
            
    print(f"The number of 3-matchings in {graph_name} is: {matching_count}")
    return matching_count

# --- Graph Definitions ---
# Both graphs are bipartite, 3-regular, on 8 vertices.

# G1: The Cube Graph (Q3)
# Vertices are numbered 0 to 7.
adj_g1 = [
    [1, 2, 4],  # 0
    [0, 3, 5],  # 1
    [0, 3, 6],  # 2
    [1, 2, 7],  # 3
    [0, 5, 6],  # 4
    [1, 4, 7],  # 5
    [2, 4, 7],  # 6
    [3, 5, 6]   # 7
]

# G2: The C8 + Chords Graph
# An 8-cycle with edges (i, i+1) plus chords connecting opposite vertices (i, i+4).
adj_g2 = [
    [1, 7, 4],  # 0
    [0, 2, 5],  # 1
    [1, 3, 6],  # 2
    [2, 4, 7],  # 3
    [3, 5, 0],  # 4
    [4, 6, 1],  # 5
    [5, 7, 2],  # 6
    [6, 0, 3]   # 7
]

# --- Calculation and Conclusion ---
print("We will test two bipartite, 3-regular graphs on 8 vertices.")
count1 = count_3_matchings("G1 (Cube Graph)", adj_g1)
count2 = count_3_matchings("G2 (C8 + Chords)", adj_g2)

print("\n--- Conclusion ---")
if count1 != count2:
    print(f"The two graphs have a different number of 3-matchings ({count1} != {count2}).")
    print("Therefore, two bipartite, d-regular graphs on n vertices do NOT necessarily have the same number of 3-matchings.")
else:
    print("The two graphs have the same number of 3-matchings.")
    print("This specific example does not disprove the statement.")

<<<No>>>