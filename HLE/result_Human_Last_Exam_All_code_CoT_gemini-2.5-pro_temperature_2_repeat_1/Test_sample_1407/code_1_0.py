def find_temporal_fixed_point_sum():
    """
    Calculates the sum of all positive integer temporal fixed points less than 100.
    A temporal fixed point k satisfies |k+k| = |k-k| within the universe's
    temporal horizon.
    """
    temporal_horizon = 48
    upper_bound = 100
    
    fixed_points = []
    
    # We need to find positive integers k < 100 such that 2*k is a multiple of the temporal_horizon.
    # The condition is (2 * k) % temporal_horizon == 0
    for k in range(1, upper_bound):
        # A fixed point k must satisfy the condition |2k| = |0| mod h(T)
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)
            
    # Calculate the sum
    total_sum = sum(fixed_points)
    
    # Format the equation string
    equation_str = " + ".join(map(str, fixed_points)) + f" = {total_sum}"
    
    print("Based on the temporal horizon h(T)=48, a temporal fixed point 'k' must be a multiple of 24.")
    print(f"The positive integer fixed points less than {upper_bound} are: {', '.join(map(str, fixed_points))}")
    print("\nThe final sum is calculated as follows:")
    print(equation_str)

find_temporal_fixed_point_sum()
<<<240>>>