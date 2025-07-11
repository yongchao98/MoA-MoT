import math

def solve_clustering_coefficient():
    """
    Calculates the clustering coefficient of the central node in a recursively
    defined graph structure.
    """

    # --- Level 1: The initial 4-node "blue triad" (a K4 graph) ---
    # In a K4, any node is connected to the other 3 nodes.
    # The central node C has 3 neighbors.
    k1 = 3
    # These 3 neighbors are all connected to each other, forming a triangle.
    # So there are 3 edges between them.
    E1 = 3
    # Properties of the base graph (K4)
    nodes_g1 = 4
    edges_g1 = 6  # A K4 graph has (4*3)/2 = 6 edges

    # --- Level 2: The 16-node structure (blue + green) ---
    # This is formed by a central G1 and 3 outer replica G1s.
    # The degree of the central node (k2) is its degree in the central G1 (k1)
    # plus connections to all nodes in the 3 outer G1s.
    k2 = k1 + 3 * nodes_g1
    # The edges between neighbors (E2) are composed of:
    # 1. Edges between neighbors from the central G1 (E1).
    # 2. Edges within each of the 3 outer G1s (3 * edges_g1).
    # 3. Edges from the 3 "base" nodes of the central G1 to the outer G1s (3 * nodes_g1).
    E2 = E1 + 3 * edges_g1 + 3 * nodes_g1
    # Properties of the G2 graph
    nodes_g2 = 4 * nodes_g1
    edges_g2 = 4 * edges_g1 + 6 * nodes_g1 # 4 internal K4s + connections between them

    # --- Level 3: The 64-node structure (blue/green + red) ---
    # This is formed by a central G2 and 3 outer replica G2s.
    # The logic is the same as for Level 2.
    k3 = k2 + 3 * nodes_g2
    E3 = E2 + 3 * edges_g2 + 3 * nodes_g2

    # --- Final Calculation for the central node of the 64-node graph ---
    k = k3
    E = E3

    # The formula for the clustering coefficient C of a node is:
    # C = (2 * E) / (k * (k - 1))
    # where k is the number of neighbors (degree) and E is the number of edges between them.
    
    numerator = 2 * E
    denominator = k * (k - 1)

    # Simplify the fraction by finding the greatest common divisor.
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    print("Calculation Steps:")
    print(f"1. For the 16-node structure, the central node has k={k2} neighbors and there are E={E2} edges between them.")
    print(f"2. For the 64-node structure, the central node has k={k} neighbors and there are E={E} edges between them.")
    print("\nFinal Clustering Coefficient Calculation:")
    print(f"Formula: C = (2 * E) / (k * (k - 1))")
    print(f"Plugging in the values: C = (2 * {E}) / ({k} * ({k} - 1))")
    print(f"C = {numerator} / {denominator}")
    print(f"Simplified Fraction: C = {simplified_numerator}/{simplified_denominator}")
    print(f"Decimal Value: C ≈ {numerator / denominator:.5f}")


solve_clustering_coefficient()
<<<25/217>>>