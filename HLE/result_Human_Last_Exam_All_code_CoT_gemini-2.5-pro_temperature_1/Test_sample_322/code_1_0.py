import math

def prove_impossibility():
    """
    This function demonstrates the mathematical impossibility of solving the
    three utilities problem on a 2D plane.
    """
    # The problem describes a complete bipartite graph K_3,3.
    # One set of vertices represents the 3 houses.
    # The other set of vertices represents the 3 utilities.
    num_houses = 3
    num_utilities = 3

    # V = total number of vertices
    V = num_houses + num_utilities

    # E = total number of edges (each house connects to each utility)
    E = num_houses * num_utilities

    print("--- Analysis of the Three Utilities Problem (Graph K_3,3) ---")
    print(f"Number of vertices (V): {num_houses} houses + {num_utilities} utilities = {V}")
    print(f"Number of edges (E): {num_houses} houses * {num_utilities} utilities = {E}")
    print("-" * 60)

    # Step 1: Assume the graph is planar and apply Euler's Formula for planar graphs.
    # Formula: V - E + F = 2, where F is the number of faces.
    print("Step 1: Assume K_3,3 is planar and use Euler's Formula (V - E + F = 2).")
    # We solve for F: F = 2 - V + E
    F_from_euler = 2 - V + E
    print(f"Solving for F: F = 2 - {V} + {E} = {F_from_euler}")
    print(f"If K_3,3 were planar, it must have exactly {F_from_euler} faces.")
    print("-" * 60)

    # Step 2: Analyze the properties of the faces in K_3,3.
    # K_3,3 is bipartite, so it has no odd-length cycles. The shortest cycle is 4.
    # Therefore, every face in a planar embedding must be bounded by at least 4 edges.
    min_edges_per_face = 4
    print("Step 2: Analyze the properties of the graph's faces.")
    print("K_3,3 is a bipartite graph, meaning it has no cycles of odd length.")
    print(f"The shortest possible cycle has 4 edges (e.g., H1-U1-H2-U2-H1).")
    print(f"Therefore, every face must be bounded by at least {min_edges_per_face} edges.")
    print("\nAn inequality relates the number of edges (E) and faces (F):")
    print("Sum of edges for all faces <= 2 * E")
    print(f"F * (minimum edges per face) <= 2 * E")
    print(f"F * {min_edges_per_face} <= 2 * {E}")

    # We solve for the maximum possible number of faces based on this property.
    F_max_from_property = (2 * E) / min_edges_per_face
    print(f"Solving for F: F <= (2 * {E}) / {min_edges_per_face} = {2 * E / min_edges_per_face}")
    print(f"So, the graph's properties require the number of faces to be at most {F_max_from_property}.")
    print("-" * 60)

    # Step 3: Check for a contradiction.
    print("Step 3: Check for a contradiction.")
    print(f"From Euler's Formula, we need F = {F_from_euler}.")
    print(f"From the graph's properties, we need F <= {F_max_from_property}.")

    is_contradiction = F_from_euler > F_max_from_property
    print(f"\nThe statement '{F_from_euler} <= {F_max_from_property}' is {not is_contradiction}.")
    if is_contradiction:
        print("\nThis is a contradiction. Our initial assumption that the graph is planar must be false.")
        print("Conclusion: It is mathematically impossible to draw K_3,3 on a 2D plane without edge crossings.")
        print("This result is formally stated by Kuratowski's Theorem.")

prove_impossibility()