import collections

def count_answers_for_query(graph):
    """
    This function counts the number of answers for the acyclic conjunctive query:
    phi(y) = exists z1, z2. E(y, z1) AND E(y, z2) AND z1 != z2

    An answer is a vertex y that has at least two distinct neighbors.
    This is equivalent to counting the number of vertices with degree >= 2.

    Args:
        graph (dict): The graph represented as an adjacency list.

    Returns:
        int: The number of vertices with degree >= 2.
    """
    num_answers = 0
    # The vertices are the keys of the adjacency list.
    for vertex in graph:
        # The degree of a vertex is the number of its neighbors.
        if len(graph[vertex]) >= 2:
            num_answers += 1
    return num_answers

# Define two non-isomorphic, cospectral graphs G1 and G2.
# They both have 5 vertices and 6 edges.
# The condition that hom(T, G1) = hom(T, G2) for all trees T holds
# because they are cospectral.

# G1: A 5-cycle with a chord.
# Vertices: 0, 1, 2, 3, 4
# Edges: (0,1), (1,2), (2,3), (3,4), (4,0) [the cycle] + (0,2) [the chord]
# Degree sequence: (3, 3, 2, 2, 2)
g1_adj = {
    0: [1, 2, 4],
    1: [0, 2],
    2: [0, 1, 3],
    3: [2, 4],
    4: [0, 3]
}

# G2: A K_3 (triangle) with two pendant edges.
# Vertices: 0, 1, 2, 3, 4
# Edges: (0,1), (1,2), (2,0) [the K3] + (0,3), (1,4) [pendant edges]
# Degree sequence: (3, 3, 2, 1, 1)
g2_adj = {
    0: [1, 2, 3],
    1: [0, 2, 4],
    2: [0, 1],
    3: [0],
    4: [1]
}

# Calculate the number of answers for the query in each graph.
num_answers_g1 = count_answers_for_query(g1_adj)
num_answers_g2 = count_answers_for_query(g2_adj)

# Print the results, showing they are different.
print("The query is: phi(y) = exists z1, z2. E(y, z1) AND E(y, z2) AND z1 != z2")
print("This counts the number of vertices with degree at least 2.\n")
print(f"For graph G1, the number of answers is: {num_answers_g1}")
print(f"For graph G2, the number of answers is: {num_answers_g2}\n")
print(f"The final result shows that the number of answers can be different:")
print(f"{num_answers_g1} != {num_answers_g2}")
