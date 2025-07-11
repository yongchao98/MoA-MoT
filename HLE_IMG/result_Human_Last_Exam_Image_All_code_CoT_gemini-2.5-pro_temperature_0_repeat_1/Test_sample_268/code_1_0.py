import math

def calculate_clustering_coefficient():
    """
    Calculates the clustering coefficient for the central node of the 64-node graph.
    """
    # The problem asks for the clustering coefficient of the central node of the 3rd level structure (64-node graph).

    # Step 1: Determine the degree 'k' of the central node.
    # The central node is the center of a 16-node (L2) substructure.
    # Its neighbors are the other 15 nodes in its own L2 substructure.
    neighbors_in_central_substructure = 15
    # It is also connected to every node in the three outer 16-node (L2) substructures.
    # Number of outer substructures = 3
    # Nodes per substructure = 16
    neighbors_in_outer_substructures = 3 * 16
    # Total degree k is the sum of these neighbors.
    k = neighbors_in_central_substructure + neighbors_in_outer_substructures

    # Step 2: Determine the number of edges 'E' between the neighbors of the central node.
    # These edges exist within the substructures.

    # Edges within the central 16-node substructure's neighbors:
    # The central L2 is made of a central L1 (K4) and 3 peripheral L1s (K4s).
    # The 3 neighbors from the central L1 form a triangle.
    edges_from_central_L1_triangle = 3
    # Each of the 3 peripheral L1s is a K4 (6 edges), and all its nodes are neighbors.
    edges_from_peripheral_L1s = 3 * 6
    edges_in_central_substructure_neighbors = edges_from_central_L1_triangle + edges_from_peripheral_L1s

    # Edges within the three outer 16-node substructures:
    # We need the total number of edges in one 16-node (L2) structure.
    # An L2 structure has 4 L1s (K4s), each with 6 edges -> 4 * 6 = 24 edges.
    # It has additional connections from the center to the 3 central nodes and 9 peripheral nodes of the outer L1s -> 12 edges.
    edges_per_L2_substructure = (4 * 6) + 12
    # Total edges from the three outer substructures.
    edges_in_outer_substructures = 3 * edges_per_L2_substructure

    # Total edges E between all neighbors.
    E = edges_in_central_substructure_neighbors + edges_in_outer_substructures

    # Step 3: Calculate the clustering coefficient C = (2 * E) / (k * (k - 1)).
    numerator = 2 * E
    denominator = k * (k - 1)

    # To simplify the fraction, find the greatest common divisor (GCD).
    common_divisor = math.gcd(numerator, denominator)
    simple_numerator = numerator // common_divisor
    simple_denominator = denominator // common_divisor

    print("Calculation of the Clustering Coefficient (C) of the Central Node")
    print("C = (2 * E) / (k * (k - 1))")
    print("-" * 30)
    print(f"Degree of the central node (k): {k}")
    print(f"Number of edges between neighbors (E): {E}")
    print("-" * 30)
    print(f"C = (2 * {E}) / ({k} * ({k} - 1))")
    print(f"C = {numerator} / ({k} * {k-1})")
    print(f"C = {numerator} / {denominator}")
    print(f"Simplified, C = {simple_numerator} / {simple_denominator}")

calculate_clustering_coefficient()
<<<43/651>>>