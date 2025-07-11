def solve_snail_problem():
    """
    Calculates the maximal distance the snail could have advanced.

    The problem is solved by decomposing the 7-minute interval into three parts:
    - A 3-minute interval from t=0 to t=3.
    - A 1-minute interval from t=3 to t=4.
    - A 3-minute interval from t=4 to t=7.

    We find the maximum possible distance for each sub-interval and sum them up.
    - L(t) is the maximum distance in a t-minute interval.
    - It can be shown that L(3) = 4 meters.
    - It can be shown that L(1) = 1 meter (by definition).

    The total distance is L(3) + L(1) + L(3).
    """
    
    # Maximum distance in the first 3-minute interval [0, 3]
    dist_part1 = 4
    
    # Maximum distance in the middle 1-minute interval [3, 4]
    dist_part2 = 1
    
    # Maximum distance in the final 3-minute interval [4, 7]
    dist_part3 = 4
    
    # Total maximal distance
    total_distance = dist_part1 + dist_part2 + dist_part3
    
    print(f"The maximal distance is the sum of the maximal distances for three sub-intervals:")
    print(f"A 3-minute part, a 1-minute part, and another 3-minute part.")
    print(f"Maximal distance = {dist_part1}m + {dist_part2}m + {dist_part3}m = {total_distance}m")

solve_snail_problem()