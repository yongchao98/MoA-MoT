def solve_temporal_fixed_points():
    """
    Finds the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point k is defined by the condition that k ⊗ k is the same
    in both forward (|k+k|) and backward (|k-k|) time-flows, considering a
    temporal horizon h(T) = 48. This implies modular arithmetic.

    The condition simplifies to (2 * k) % 48 == 0.
    """
    temporal_horizon = 48
    limit = 100
    fixed_points = []

    # Find all positive integer fixed points less than the limit
    for k in range(1, limit):
        # The condition for a fixed point is |2k| = |0| within the temporal horizon.
        # This translates to (2 * k) % temporal_horizon == 0.
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)

    # Calculate the sum
    total_sum = sum(fixed_points)

    # Format the output string to show the equation
    equation_str = " + ".join(map(str, fixed_points)) + f" = {total_sum}"

    print("The positive integer temporal fixed points less than 100 are the numbers k where (2 * k) is a multiple of 48.")
    print("Found fixed points and their sum:")
    print(equation_str)

solve_temporal_fixed_points()