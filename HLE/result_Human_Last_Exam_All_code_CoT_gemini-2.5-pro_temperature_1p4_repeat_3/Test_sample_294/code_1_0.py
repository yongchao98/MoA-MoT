import math

def combinations(n, k):
    """Calculates n choose k."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def count_2_matchings(graph):
    """
    Calculates the number of 2-matchings in a graph.
    The graph is represented as an adjacency list (dict of lists).
    The formula used is: (# 2-matchings) = (m choose 2) - sum(deg(v) choose 2 for v in V)
    where m is the number of edges and deg(v) is the degree of vertex v.
    """
    if not graph:
        print("Graph is empty.")
        return

    # Calculate the number of vertices (n) and edges (m)
    num_vertices = len(graph)
    degrees = {v: len(neighbors) for v, neighbors in graph.items()}
    # Each edge is counted twice in the sum of degrees
    num_edges = sum(degrees.values()) // 2

    print(f"Graph has {num_vertices} vertices and {num_edges} edges.")

    # Term 1: Total number of pairs of edges
    total_edge_pairs = combinations(num_edges, 2)
    print("\nEquation: (# 2-matchings) = (m choose 2) - sum(deg(v) choose 2)")
    print("---------------------------------------------------------------")
    print(f"Number of edges m = {num_edges}")
    print(f"Term 1 (m choose 2) = combinations({num_edges}, 2) = {total_edge_pairs}")

    # Term 2: Sum of pairs of incident edges
    incident_edge_pairs_sum = 0
    calculation_steps = []
    for vertex, degree in degrees.items():
        inc_pairs = combinations(degree, 2)
        incident_edge_pairs_sum += inc_pairs
        if inc_pairs > 0:
            calculation_steps.append(f"deg({vertex})={degree} -> ({degree} choose 2)={inc_pairs}")

    print(f"\nTerm 2 (sum over vertices v of (deg(v) choose 2)):")
    if not calculation_steps:
        print("  All vertices have degree < 2. Sum is 0.")
    else:
        print("  " + "\n  ".join(calculation_steps))
    print(f"  Total sum = {incident_edge_pairs_sum}")

    # Final result
    num_2_matchings = total_edge_pairs - incident_edge_pairs_sum

    print("\nFinal Calculation:")
    print(f"  {total_edge_pairs} (total edge pairs)")
    print(f"- {incident_edge_pairs_sum} (incident edge pairs)")
    print("--------------------")
    print(f"  = {num_2_matchings} (2-matchings)")
    return num_2_matchings


if __name__ == "__main__":
    # Example: A cycle graph with 5 vertices (C5)
    # The 2-matchings are: {(0,1),(2,3)}, {(0,1),(3,4)}, {(1,2),(3,4)}, {(1,2),(4,0)}, {(2,3),(4,0)}
    # There are 5 such 2-matchings.
    c5_graph = {
        0: [1, 4],
        1: [0, 2],
        2: [1, 3],
        3: [2, 4],
        4: [3, 0]
    }
    print("### Counting 2-matchings in a 5-cycle graph (C5) ###")
    count_2_matchings(c5_graph)