def calculate_ghtw_for_cyclic_hypergraph(K):
    """
    Calculates the generalised hypertreewidth (g-htw) for a constructed
    3-edge hypergraph to demonstrate its unbounded nature.

    The hypergraph is constructed from three disjoint sets of vertices A, B, and C.
    We assume |A| = |B| = |C| = K.
    The hyperedges are e1 = A U B, e2 = B U C, e3 = C U A.
    The g-htw of this hypergraph is |A| + |B| + |C| - 1.

    Args:
        K (int): The size for each of the disjoint vertex sets A, B, and C.
                 A larger K demonstrates how the g-htw can grow without bound.
    """
    if not isinstance(K, int) or K <= 0:
        print("Error: K must be a positive integer.")
        return

    size_A = K
    size_B = K
    size_C = K

    # The g-htw for this construction is the sum of the sizes of the
    # disjoint sets minus one.
    result = size_A + size_B + size_C - 1

    print("For a constructed hypergraph with 3 hyperedges and component sets of size K:")
    print(f"K = {K}")
    print(f"The size of set A, |A| = {size_A}")
    print(f"The size of set B, |B| = {size_B}")
    print(f"The size of set C, |C| = {size_C}")
    print("\nThe g-htw is calculated using the formula: |A| + |B| + |C| - 1")
    # Final output showing the numbers in the equation
    print(f"Final Equation: {size_A} + {size_B} + {size_C} - 1 = {result}")
    print(f"\nBy choosing a larger K, we can make the g-htw arbitrarily large.")

# Demonstrate with an example value for K.
# The user can change this value to see the g-htw change.
example_K = 100
calculate_ghtw_for_cyclic_hypergraph(example_K)
