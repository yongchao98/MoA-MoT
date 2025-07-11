def solve_temporal_puzzle():
    """
    Finds the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point k satisfies the condition:
    |k + k| is equivalent to |k - k| within the universal temporal horizon.
    The temporal horizon h(T) is 48.
    This translates to the mathematical condition: (2 * k) % 48 == 0.
    """
    h_T = 48
    limit = 100
    fixed_points = []

    # Iterate through all positive integers k less than the limit
    for k in range(1, limit):
        # Condition for a temporal fixed point: 2*k must be a multiple of h_T
        if (2 * k) % h_T == 0:
            fixed_points.append(k)

    # Calculate the sum
    total_sum = sum(fixed_points)

    # Format the equation string as requested
    equation_str = " + ".join(map(str, fixed_points)) + f" = {total_sum}"

    print("The positive integer temporal fixed points less than 100 are:")
    print(fixed_points)
    print("\nThe sum of these fixed points is calculated as follows:")
    print(equation_str)


solve_temporal_puzzle()
<<<240>>>