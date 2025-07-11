def solve_temporal_fixed_points():
    """
    Solves for the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point k is a value where the result of the temporal
    multiplication k ⊗ k is the same regardless of time-flow direction.
    - Forward flow: k ⊗ k = |k + k| = |2k|
    - Backward flow: k ⊗ k = |k - k| = |0| = 0

    The "temporal horizon" h(T) = 48 implies all calculations are modulo 48.
    So, the condition for a fixed point k is: |2k| mod 48 = 0.
    For positive k, this simplifies to (2 * k) % 48 == 0.
    This means 2k must be a multiple of 48, so k must be a multiple of 24.

    The program will find all multiples of 24 that are positive and less than 100,
    and then compute their sum.
    """

    temporal_horizon = 48
    limit = 100
    fixed_points = []

    # Iterate through all positive integers less than the limit
    for k in range(1, limit):
        # Check if k satisfies the temporal fixed point condition
        # (2 * k) must be a multiple of the temporal_horizon
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)

    # Calculate the sum of the found fixed points
    total_sum = sum(fixed_points)

    # Format the equation string showing each number
    equation_str = " + ".join(map(str, fixed_points))

    # Print the final result including the full equation
    print("The sum of all positive integer temporal fixed points less than 100 is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>