def solve_snail_problem():
    """
    Calculates the maximal distance the snail could have advanced in 7 minutes.
    """
    L = 7  # Total duration in minutes

    # The problem is a known mathematical puzzle. The maximal distance D for a
    # duration of L minutes follows the formula:
    # D = 3 * floor(L / 2) + (L mod 2)
    # This is derived by breaking the L-minute interval into floor(L/2)
    # 2-minute segments and one 1-minute segment if L is odd.
    # In each 2-minute segment, a maximum of 3 meters can be traveled.
    # In the remaining 1-minute segment, a maximum of 1 meter can be traveled.

    num_2_min_segments = L // 2
    num_1_min_segments = L % 2

    dist_from_2_min = num_2_min_segments * 3
    dist_from_1_min = num_1_min_segments * 1

    total_distance = dist_from_2_min + dist_from_1_min

    print("The total duration of observation is {} minutes.".format(L))
    print("This can be broken down into {} segments of 2 minutes and {} segment of 1 minute.".format(num_2_min_segments, num_1_min_segments))
    print("Maximal distance in a 2-minute segment: 3 meters")
    print("Maximal distance in a 1-minute segment: 1 meter")
    print("The total maximal distance is ({} * 3) + ({} * 1) = {} meters.".format(num_2_min_segments, num_1_min_segments, total_distance))

solve_snail_problem()