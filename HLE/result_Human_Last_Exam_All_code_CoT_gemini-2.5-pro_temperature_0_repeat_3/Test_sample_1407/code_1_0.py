def solve_temporal_fixed_points():
    """
    Finds the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point k satisfies the condition |k + k| = |k - k|
    within the context of a temporal horizon h(T) = 48. This implies
    that the operations are performed modulo 48.

    The condition simplifies to (2 * k) % 48 == 0.
    """
    temporal_horizon = 48
    limit = 100
    fixed_points = []

    # We are looking for positive integers k < 100
    for k in range(1, limit):
        # The condition for a fixed point is |2k| mod 48 = |0| mod 48
        # Since k is positive, this is (2 * k) % 48 == 0
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)

    # Calculate the sum
    total_sum = sum(fixed_points)

    # Create the equation string
    equation_str = " + ".join(map(str, fixed_points))

    # Print the results
    print(f"The positive integer temporal fixed points less than {limit} are: {fixed_points}")
    print("The sum is calculated as:")
    print(f"{equation_str} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>