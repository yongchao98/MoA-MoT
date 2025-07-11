import math

def count_2_matchings(graph):
    """
    Calculates the number of 2-matchings in a graph.

    A 2-matching is a set of 2 edges that do not share any vertices.
    The number is calculated using the principle of inclusion-exclusion:
    Total number of pairs of edges - Number of pairs of edges that share a vertex.

    Args:
        graph (dict): A dictionary representing the graph as an adjacency list.
                      Keys are vertices (int), values are lists of neighbors.
    """
    num_vertices = len(graph)
    if num_vertices == 0:
        print("Graph is empty.")
        return

    # Step 1: Calculate the total number of edges, |E|
    num_edges = 0
    degrees = {}
    for vertex, neighbors in graph.items():
        degree = len(neighbors)
        degrees[vertex] = degree
        num_edges += degree
    # Each edge is counted twice (once for each endpoint), so divide by 2
    num_edges //= 2

    # Step 2: Calculate the total number of pairs of edges, C(|E|, 2)
    total_edge_pairs = math.comb(num_edges, 2)

    # Step 3: Calculate the number of conflicting edge pairs.
    # A conflicting pair of edges shares a vertex, forming a "cherry" or path of length 2.
    # The number of such paths centered at a vertex v is C(degree(v), 2).
    num_cherries = 0
    for vertex in graph:
        degree = degrees[vertex]
        if degree >= 2:
            num_cherries += math.comb(degree, 2)
            
    # Step 4: Calculate the number of 2-matchings
    num_2_matchings = total_edge_pairs - num_cherries

    # Output the final equation with the computed numbers
    print(f"Graph details:")
    print(f"  - Number of vertices: {num_vertices}")
    print(f"  - Number of edges |E|: {num_edges}")
    
    print("\nCalculation for the number of 2-matchings:")
    print(f"  - Total pairs of edges = C(|E|, 2) = C({num_edges}, 2) = {total_edge_pairs}")
    print(f"  - Pairs of edges sharing a vertex (cherries) = {num_cherries}")
    print("\nFinal Equation:")
    print(f"Number of 2-matchings = C({num_edges}, 2) - {num_cherries} = {total_edge_pairs} - {num_cherries} = {num_2_matchings}")


if __name__ == '__main__':
    # Example Graph: A pentagon (C5) and a disconnected edge (K2)
    # Vertices 0-4 form the C5, vertices 5-6 form the K2.
    # The 2-matchings in C5 are: {(0,1),(2,3)}, {(1,2),(3,4)}, {(2,3),(4,0)}, {(3,4),(0,1)}, {(4,0),(1,2)} -> 5
    # Additional 2-matchings involve the edge (5,6) paired with any of the 5 edges of C5 -> 5
    # Total = 5 + 5 = 10
    
    example_graph = {
        0: [1, 4],
        1: [0, 2],
        2: [1, 3],
        3: [2, 4],
        4: [0, 3],
        5: [6],
        6: [5]
    }
    count_2_matchings(example_graph)
