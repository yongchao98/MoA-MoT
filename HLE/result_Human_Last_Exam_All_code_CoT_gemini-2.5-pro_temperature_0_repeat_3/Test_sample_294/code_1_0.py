import math

def count_2_matchings(adj_list):
    """
    Counts the number of 2-matchings in a graph represented by an adjacency list.
    A 2-matching is a set of 2 edges that do not share any vertices.

    The algorithm runs in O(|V|^2) time, which is subcubic.

    The formula used is:
    Number of 2-matchings = (Total pairs of edges) - (Pairs of edges sharing a vertex)
                          = (|E| choose 2) - sum_{v in V} (degree(v) choose 2)
    """
    num_vertices = len(adj_list)
    if num_vertices == 0:
        print("Graph is empty.")
        return

    degrees = [len(adj_list[i]) for i in range(num_vertices)]
    num_edges = sum(degrees) // 2

    print(f"Graph has {num_vertices} vertices and {num_edges} edges.")
    print("-" * 20)
    print("Equation: N(2-matchings) = (|E| C 2) - SUM(d(v) C 2)")
    print("-" * 20)

    # Helper function for combinations "n choose 2"
    def n_choose_2(n):
        if n < 2:
            return 0
        return n * (n - 1) // 2

    # Calculate the first term: |E| choose 2
    total_edge_pairs = n_choose_2(num_edges)
    print(f"Term 1: Total number of edge pairs, |E| C 2")
    print(f"|E| = {num_edges}")
    print(f"|E| C 2 = {total_edge_pairs}")
    print("-" * 20)

    # Calculate the second term: sum over all vertices v of (degree(v) choose 2)
    # This counts pairs of edges that are incident to the same vertex.
    colliding_edge_pairs = sum(n_choose_2(d) for d in degrees)
    print(f"Term 2: Number of colliding edge pairs, SUM(d(v) C 2)")
    print(f"Degrees d(v): {degrees}")
    print(f"SUM(d(v) C 2) = {colliding_edge_pairs}")
    print("-" * 20)

    # Final calculation
    num_2_matchings = total_edge_pairs - colliding_edge_pairs

    print("Final Calculation:")
    print(f"Number of 2-matchings = {total_edge_pairs} - {colliding_edge_pairs}")
    print(f"Result = {num_2_matchings}")


if __name__ == '__main__':
    # Example Graph: A complete graph K4
    # Vertices: 0, 1, 2, 3
    # Edges: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    # The 2-matchings are: {(0,1),(2,3)}, {(0,2),(1,3)}, {(0,3),(1,2)}. Total = 3.
    k4_adj_list = [
        [1, 2, 3],  # Neighbors of vertex 0
        [0, 2, 3],  # Neighbors of vertex 1
        [0, 1, 3],  # Neighbors of vertex 2
        [0, 1, 2]   # Neighbors of vertex 3
    ]
    print("--- Counting 2-matchings in a complete graph K4 ---")
    count_2_matchings(k4_adj_list)
