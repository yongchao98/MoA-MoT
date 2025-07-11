def solve_temporal_fixed_points():
    """
    Finds the sum of all positive integer temporal fixed points less than 100.

    A temporal fixed point 'k' is defined by the condition that its temporal multiplication
    yields the same result in both forward and backward time-flows.
    - Forward: |k + k|
    - Backward: |k - k|

    The universe has a "temporal horizon" h(T) = 48, which we interpret as
    a modulus for all temporal calculations.

    The condition becomes: |2k| mod 48 == |0| mod 48
    For positive integers k, this simplifies to: (2 * k) % 48 == 0.
    This means k must be a multiple of 24.
    """
    
    limit = 100
    temporal_horizon = 48
    
    fixed_points = []
    
    # We are looking for k such that (2*k) is a multiple of 48.
    # This is equivalent to k being a multiple of 24.
    # We can iterate through multiples of 24 directly.
    k = 24
    while k < limit:
        fixed_points.append(k)
        k += 24
        
    total_sum = sum(fixed_points)
    
    # Create the equation string for printing
    equation_str = " + ".join(map(str, fixed_points))
    
    print(f"{equation_str} = {total_sum}")

solve_temporal_fixed_points()