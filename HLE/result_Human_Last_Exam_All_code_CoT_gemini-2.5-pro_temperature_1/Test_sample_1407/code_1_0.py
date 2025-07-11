def solve_temporal_fixed_points():
    """
    Finds and sums all positive integer temporal fixed points less than 100.
    """
    # The temporal horizon, h(T), is a universal constant acting as a modulus.
    temporal_horizon = 48

    # A temporal fixed point 'k' is where |(k+k)| = |(k-k)|.
    # With the temporal horizon, this becomes (2*k) % h(T) == 0.
    
    fixed_points = []
    # We are looking for positive integers less than 100.
    for k in range(1, 100):
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)
            
    # Calculate the sum of the found fixed points.
    total_sum = sum(fixed_points)

    # Format the numbers into an equation string.
    equation_str = " + ".join(map(str, fixed_points))
    
    print(f"The positive integer temporal fixed points less than 100 are: {fixed_points}")
    print("The sum is found by solving the equation:")
    print(f"{equation_str} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>