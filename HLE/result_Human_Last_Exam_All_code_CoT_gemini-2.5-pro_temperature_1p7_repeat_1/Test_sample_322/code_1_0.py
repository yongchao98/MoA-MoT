def check_k33_planarity():
    """
    Demonstrates the impossibility of solving the three utilities problem
    by using Euler's formula for planar graphs to prove K3,3 is non-planar.
    """

    # 1. Define the properties of the K3,3 graph.
    # 3 houses + 3 utilities
    V = 6
    # 3 houses * 3 utilities
    E = 9

    print("--- Proof of Impossibility for the Three Utilities Problem (K3,3 Graph) ---")
    print(f"Number of Vertices (V): {V}")
    print(f"Number of Edges (E): {E}")
    print("-" * 65)

    # 2. Use Euler's formula for planar graphs: V - E + F = 2
    # If K3,3 were planar, we can calculate the number of faces (F) it would have.
    # F = E - V + 2
    F = E - V + 2

    print("Step 1: Applying Euler's Formula (V - E + F = 2)")
    print("If the graph were planar, the number of faces (F) it would create is:")
    print(f"F = {E} - {V} + 2 = {F}")
    print("-" * 65)

    # 3. Use the property that K3,3 is bipartite.
    # This means it has no odd-length cycles. The shortest possible cycle has 4 edges.
    # Therefore, every face in a planar drawing must be bounded by at least 4 edges.
    min_cycle_length = 4
    print("Step 2: Using the Bipartite Graph Property")
    print(f"K3,3 is a bipartite graph, so the shortest cycle length is {min_cycle_length}.")
    print("This means every face must be bounded by at least 4 edges.")
    print("-" * 65)

    # 4. Check the condition for planar graphs.
    # For any planar graph, the sum of the lengths of all face boundaries must equal 2*E.
    # This leads to the inequality: min_cycle_length * F <= 2 * E
    # We will now check if this holds true for K3,3.

    lhs = min_cycle_length * F
    rhs = 2 * E

    print("Step 3: Checking the Planarity Condition (min_cycle_length * F <= 2 * E)")
    print(f"We must check if: {min_cycle_length} * {F} <= 2 * {E}")
    print(f"Calculating each side: {lhs} <= {rhs}")

    print("\n--- Conclusion ---")
    if lhs <= rhs:
        print(f"The inequality {lhs} <= {rhs} is true. This proof is inconclusive.")
    else:
        print(f"The inequality {lhs} <= {rhs} is FALSE.")
        print("This is a contradiction. The assumption that K3,3 is planar must be incorrect.")
        print("\nTherefore, it is mathematically impossible to connect all three houses to all three utilities on a 2D plane without any lines crossing.")

# Run the proof
check_k33_planarity()
