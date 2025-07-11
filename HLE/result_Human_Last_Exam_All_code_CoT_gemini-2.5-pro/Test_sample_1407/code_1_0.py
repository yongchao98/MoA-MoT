def find_temporal_fixed_point_sum():
    """
    Finds the sum of all positive integer temporal fixed points less than 100.
    
    A fixed point k satisfies the condition |k+k| = |k-k| within the universe's
    arithmetic, which is governed by the temporal horizon h(T) = 48 acting as a modulus.
    The condition is therefore (2 * k) % 48 == 0.
    """
    
    # The temporal horizon constant
    h = 48
    
    # The upper limit for our search (less than 100)
    limit = 100
    
    # List to store the found fixed points
    fixed_points = []
    
    # Iterate through all positive integers less than the limit
    for k in range(1, limit):
        # A fixed point k must satisfy (2*k) mod h = 0
        if (2 * k) % h == 0:
            fixed_points.append(k)
            
    # Calculate the sum of the found fixed points
    total_sum = sum(fixed_points)
    
    print(f"The positive integer temporal fixed points less than {limit} are: {', '.join(map(str, fixed_points))}")
    
    # Create the equation string as required
    equation_str = " + ".join(map(str, fixed_points))
    
    print(f"The sum is calculated as: {equation_str} = {total_sum}")

find_temporal_fixed_point_sum()
<<<240>>>