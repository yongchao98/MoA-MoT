import itertools

def count_3_matchings(edges):
    """
    Counts the number of 3-matchings in a graph.
    A 3-matching is a set of 3 edges with no common vertices.
    """
    count = 0
    # Iterate through all combinations of 3 edges
    for edge_triplet in itertools.combinations(edges, 3):
        # A 3-matching consists of 3 edges and 6 distinct vertices
        vertices = set()
        for edge in edge_triplet:
            vertices.add(edge[0])
            vertices.add(edge[1])
        
        # If there are 6 unique vertices, the 3 edges are disjoint
        if len(vertices) == 6:
            count += 1
    return count

# Graph G1: Two disjoint copies of K_{3,3}
# Vertices 0-5 form the first K_{3,3}. Partitions: {0,1,2} and {3,4,5}
# Vertices 6-11 form the second K_{3,3}. Partitions: {6,7,8} and {9,10,11}
g1_edges = []
# First K_{3,3}
for i in range(3):
    for j in range(3, 6):
        g1_edges.append(tuple(sorted((i, j))))
# Second K_{3,3}
for i in range(6, 9):
    for j in range(9, 12):
        g1_edges.append(tuple(sorted((i, j))))

# Graph G2: A connected, 3-regular, bipartite graph on 12 vertices
# Partitions: U = {0,...,5} and V = {6,...,11}
# Edge definition: vertex i in U is connected to (i+6), ((i-1)%6 + 6), and ((i+2)%6 + 6) in V
g2_edges = []
for i in range(6):
    # Edges from partition U to V
    v1 = i + 6
    v2 = (i - 1 + 6) % 6 + 6
    v3 = (i + 2) % 6 + 6
    g2_edges.append(tuple(sorted((i, v1))))
    g2_edges.append(tuple(sorted((i, v2))))
    g2_edges.append(tuple(sorted((i, v3))))

# Calculate the number of 3-matchings for each graph
num_matchings_g1 = count_3_matchings(g1_edges)
num_matchings_g2 = count_3_matchings(g2_edges)

print("Both graphs are bipartite, 3-regular, and have 12 vertices.")
print(f"Graph G1 (2 disjoint K_3,3) has {num_matchings_g1} 3-matchings.")
print(f"Graph G2 (a connected graph) has {num_matchings_g2} 3-matchings.")

if num_matchings_g1 != num_matchings_g2:
    print("\nSince the number of 3-matchings is different, the two graphs do not necessarily have the same number of 3-matchings.")
else:
    print("\nThe number of 3-matchings is the same for this example.")
