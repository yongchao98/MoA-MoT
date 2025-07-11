def solve_temporal_fixed_points():
    """
    Finds and sums all positive integer temporal fixed points less than 100.
    """
    # Define the constants from the problem.
    temporal_horizon = 48
    limit = 100

    fixed_points = []
    
    # A temporal fixed point 'k' satisfies the condition:
    # |k + k| is equivalent to |k - k| within the universal rules.
    # Forward result: |2k| = 2k (for positive k)
    # Backward result: |k - k| = 0
    # The "temporal horizon" h(T) implies a modulus operation.
    # The condition becomes: 2k is congruent to 0, modulo the temporal horizon.
    # So, (2 * k) % temporal_horizon == 0
    
    # Iterate through all positive integers less than the limit.
    for k in range(1, limit):
        # Check if 2*k is a multiple of the temporal horizon.
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)
            
    # Calculate the sum of all found fixed points.
    total_sum = sum(fixed_points)
    
    # Format the output string to show the equation.
    equation_str = " + ".join(map(str, fixed_points))
    
    # Print the final equation with the sum.
    print(f"The temporal fixed points less than {limit} are the positive integers k such that (2 * k) % {temporal_horizon} == 0.")
    print("These points are: ", ", ".join(map(str, fixed_points)))
    print("\nCalculating the sum:")
    print(f"{equation_str} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>