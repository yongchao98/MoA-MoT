def find_temporal_fixed_point_sum():
    """
    Finds all positive integer temporal fixed points less than 100 and calculates their sum.

    A temporal fixed point 'k' is defined by the condition that the results of its
    "temporal multiplication" are the same regardless of time-flow direction, considering
    a "temporal horizon" h(T) = 48.

    Forward operation: |k + k|
    Backward operation: |k - k|
    Temporal Horizon, h(T): 48 (implies calculations are modulo 48)

    The condition for a fixed point 'k' is: |k + k| mod 48 == |k - k| mod 48
    For a positive integer k, this simplifies to: 2*k % 48 == 0.
    """
    
    temporal_horizon = 48
    limit = 100
    fixed_points = []

    # Iterate through all positive integers less than the limit
    for k in range(1, limit):
        # A temporal fixed point k satisfies the condition: 2*k is a multiple of the temporal horizon.
        if (2 * k) % temporal_horizon == 0:
            fixed_points.append(k)
    
    # Calculate the sum of the found fixed points
    total_sum = sum(fixed_points)

    # Format the numbers into a string equation "n1 + n2 + ... = sum"
    equation_str = " + ".join(map(str, fixed_points))
    
    print(f"The temporal fixed points less than {limit} are: {', '.join(map(str, fixed_points))}")
    print("The sum is calculated as:")
    print(f"{equation_str} = {total_sum}")

find_temporal_fixed_point_sum()
<<<240>>>