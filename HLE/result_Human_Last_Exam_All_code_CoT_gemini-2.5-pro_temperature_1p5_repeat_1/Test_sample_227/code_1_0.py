def solve_snail_puzzle():
    """
    This function calculates the maximal distance the snail could have traveled.
    The problem is solved by partitioning the 7-minute interval and maximizing
    the distance in each part.

    1.  Intervals [0, 3] and [4, 7]:
        Observers watching [0,1], [1,2], [2,3] force the snail to travel 3 meters.
        p(1)-p(0)=1, p(2)-p(1)=1, p(3)-p(2)=1 => p(3)-p(0)=3.
        Similarly, observers on [4,5], [5,6], [6,7] force 3 meters of travel.
        p(5)-p(4)=1, p(6)-p(5)=1, p(7)-p(6)=1 => p(7)-p(4)=3.

    2.  Interval [3, 4]:
        The distance traveled in the middle interval [3, 4] can be maximized.
        By using a set of 4 observers that straddle the [3,4] interval
        (e.g., [2.5, 3.5], [2.75, 3.75], [3, 4], [3.25, 4.25]), it's possible for the snail
        to travel 4 meters in this 1-minute slot while each observer still
        only records 1 meter of total advancement.

    3.  Total Distance:
        The total distance is the sum of distances from these three parts.
    """
    dist_0_to_3 = 3
    dist_4_to_7 = 3
    dist_3_to_4 = 4 # The maximum possible distance in the middle interval.

    total_distance = dist_0_to_3 + dist_3_to_4 + dist_4_to_7
    
    print(f"The total distance is the sum of the distances in three segments:")
    print(f"Segment 1 (0 to 3 mins): {dist_0_to_3} meters")
    print(f"Segment 2 (3 to 4 mins): {dist_3_to_4} meters")
    print(f"Segment 3 (4 to 7 mins): {dist_4_to_7} meters")
    print(f"The equation for the final answer is: {dist_0_to_3} + {dist_3_to_4} + {dist_4_to_7} = {total_distance}")
    print(f"Maximal distance: {total_distance} meters")

solve_snail_puzzle()