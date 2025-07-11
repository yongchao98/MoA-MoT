def solve_k_vector_problem():
    """
    This script determines the smallest k for which any bridgeless 3-regular graph
    admits a valid k-vector based on established graph theory results.
    """

    # The problem defines a valid k-vector 'x' for a graph G. This vector must satisfy two conditions:
    # 1. B*x = 0, where B is the incidence matrix. For a 3-regular graph, this means at any
    #    vertex with incident edges e1, e2, e3, the assigned values must sum to zero:
    #    x_1 + x_2 + x_3 = 0.
    # 2. The entries x_i belong to the set {+/-1, +/-2, ..., +/-(k-1)}.
    # This is known as a nowhere-zero k-flow. We need the smallest k that works for *any*
    # bridgeless 3-regular graph, regardless of its specific structure (e.g., having 20 vertices).
    # We test values of k, starting from the smallest possible.

    # --- Check k = 2 ---
    # Vector entries are in {+1, -1}.
    # The equation x_1 + x_2 + x_3 = 0 cannot be satisfied with three numbers from {+1, -1}.
    # For example, 1+1-1=1 != 0. So, k=2 is not sufficient.
    k_is_sufficient = False

    # --- Check k = 3 ---
    # Vector entries are in {+/-1, +/-2}. This is a 3-flow.
    # For a 3-regular graph, having a 3-flow is equivalent to being 3-edge-colorable.
    # However, there are bridgeless 3-regular graphs (called 'snarks', e.g., the Petersen graph)
    # that are not 3-edge-colorable. Thus, k=3 is not sufficient for all graphs in the class.
    k_is_sufficient = False

    # --- Check k = 4 ---
    # Vector entries are in {+/-1, +/-2, +/-3}. This is a 4-flow.
    # The Petersen graph is a well-known counterexample. It is a bridgeless 3-regular graph,
    # and it is proven that it does not have a 4-flow.
    # Therefore, k=4 is not sufficient for all graphs in the class.
    k_is_sufficient = False

    # --- Check k = 5 ---
    # Vector entries are in {+/-1, +/-2, +/-3, +/-4}. This is a 5-flow.
    # It is a proven theorem that every bridgeless 3-regular graph has a 5-flow.
    # This result holds for all graphs in the specified class.
    # Therefore, k=5 is the smallest sufficient value.
    final_k = 5
    k_is_sufficient = True

    print(f"The smallest value of k is {final_k}.")

    # To fulfill the prompt's request, we show the numbers in the final equation.
    # For k=5, a valid flow can be constructed. For any given vertex, the equation
    # x_1 + x_2 + x_3 = 0 must hold, where x_i are from {+/-1, +/-2, +/-3, +/-4}.
    # A valid assignment at one vertex could be:
    val1 = 2
    val2 = 2
    val3 = -4
    
    print("\nThe equation at each vertex must be satisfied. For example:")
    print(f"{val1} + {val2} + ({val3}) = {val1 + val2 + val3}")

solve_k_vector_problem()