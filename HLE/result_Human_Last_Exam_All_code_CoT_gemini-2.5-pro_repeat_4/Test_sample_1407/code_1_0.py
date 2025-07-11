def solve_temporal_fixed_points():
    """
    This function finds all positive integer temporal fixed points less than a limit,
    based on the given definitions, and calculates their sum.

    A temporal fixed point 'k' is defined by the condition:
    |k + k| mod h(T) = |k - k| mod h(T)

    Given:
    - h(T) = 48 (the temporal horizon, used as a modulus)
    - The limit for k is 100.

    The condition simplifies to:
    |2k| mod 48 = |0| mod 48
    2k mod 48 = 0 (since k must be a positive integer)

    This means 2k must be a multiple of 48, which implies k must be a
    multiple of 24.
    """
    temporal_horizon = 48
    limit = 100
    
    # We need to find multiples of (temporal_horizon / 2)
    base_multiple = temporal_horizon // 2
    
    fixed_points = []
    # Iterate through positive integers less than the limit
    for k in range(1, limit):
        # A fixed point k must be a multiple of 24
        if k % base_multiple == 0:
            fixed_points.append(k)
            
    # Calculate the sum of the found fixed points
    total_sum = sum(fixed_points)
    
    # Format the numbers into an equation string, as requested
    equation_string = " + ".join(map(str, fixed_points))
    
    print(f"The positive integer temporal fixed points less than {limit} are the multiples of {base_multiple}.")
    print(f"Found fixed points: {fixed_points}")
    print("The sum is calculated as:")
    print(f"{equation_string} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>