def solve_utility_problem():
    """
    Analyzes the Three Utilities Problem using graph theory principles
    to demonstrate its impossibility.
    """
    # Step 1: Define the properties of the graph (K_3,3).
    # V = number of vertices (3 houses + 3 utilities)
    # E = number of edges (each of the 3 houses connects to each of the 3 utilities)
    V = 3 + 3
    E = 3 * 3

    print("--- Analysis of the Three Utilities Problem (Graph K_3,3) ---")
    print(f"Number of Vertices (V): {V}")
    print(f"Number of Edges (E): {E}\n")

    # Step 2: Assume the graph is planar and use Euler's formula (V - E + F = 2)
    # to find the required number of faces (F).
    print("--- Step 1: Applying Euler's Formula for Planar Graphs ---")
    print("Formula: V - E + F = 2")
    print("If the graph were planar, we can solve for the number of faces (F):")
    print(f"F = 2 - V + E")
    # Calculate F
    F = 2 - V + E
    print(f"F = 2 - {V} + {E} = {F}\n")

    # Step 3: Analyze the properties of the faces.
    # The graph is bipartite, so the shortest cycle is length 4.
    # Therefore, every face must be bounded by at least 4 edges.
    min_edges_per_face = 4
    print("--- Step 2: Analyzing the Faces of the Graph ---")
    print(f"The graph is bipartite, so the shortest cycle length is {min_edges_per_face}.")
    print("This implies every face in a planar drawing must be bounded by at least 4 edges.\n")


    # Step 4: Use this property to form an inequality and check for a contradiction.
    # For any planar graph, 2 * E must be greater than or equal to F * min_edges_per_face.
    print("--- Step 3: Checking for a Mathematical Contradiction ---")
    print("This leads to the inequality: (F * min_edges_per_face) <= (2 * E)")
    
    # Calculate both sides of the inequality
    left_side = F * min_edges_per_face
    right_side = 2 * E

    print(f"Substituting the values: ({F} * {min_edges_per_face}) <= (2 * {E})")
    print(f"Which simplifies to: {left_side} <= {right_side}")
    
    # Conclusion
    print("\n--- Conclusion ---")
    if left_side <= right_side:
        print("The inequality holds, which does not lead to a contradiction.")
    else:
        print(f"The result {left_side} <= {right_side} is FALSE.")
        print("This contradiction proves that the initial assumption was incorrect.")
        print("The graph K_3,3 is non-planar, and the problem is impossible to solve under the given constraints.")
        print("\nThe correct option is E, which states the problem is impossible due to Kuratowski's Theorem.")

solve_utility_problem()