def solve_temporal_fixed_points():
    """
    Calculates the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point k satisfies the condition:
    |k + k| mod h(T) = |k - k| mod h(T)
    where the temporal horizon h(T) = 48.

    This simplifies to (2 * k) % 48 == 0, meaning k must be a multiple of 24.
    """
    temporal_horizon = 48
    limit = 100
    fixed_points = []

    # k must be a multiple of 24.
    base_k = temporal_horizon // 2  # 48 / 2 = 24

    # Find all multiples of 24 that are less than the limit.
    current_multiple = base_k
    while current_multiple < limit:
        if current_multiple > 0:
            fixed_points.append(current_multiple)
        current_multiple += base_k

    # Calculate the sum of the found fixed points.
    total_sum = sum(fixed_points)

    # Create the equation string for display.
    # The requirement is to show each number in the final equation.
    equation_str = " + ".join(map(str, fixed_points)) + f" = {total_sum}"

    print(f"The positive integer temporal fixed points less than {limit} are: {fixed_points}")
    print("The sum is calculated as follows:")
    print(equation_str)

    # The final answer format as requested
    print(f"<<<{total_sum}>>>")

solve_temporal_fixed_points()