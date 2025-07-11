def solve_snail_problem():
    """
    Calculates the maximal distance the snail could have advanced.

    The problem can be solved by splitting the 7-minute interval into two
    overlapping sub-intervals and finding the maximal distance for each.

    Let D(T) be the maximal distance for a T-minute interval.
    For T >= 3, it can be shown that D(T) = T + 1.
    For T = 2, D(2) = 2.
    For T = 1, D(1) = 1.

    We split the 7-minute interval into two overlapping intervals:
    1. A 3-minute interval from t=0 to t=3.
    2. A 5-minute interval from t=2 to t=7.
    """

    # Let d(t) be the distance at time t. We set d(0) = 0.
    d_0 = 0

    # For the first interval [0, 3], the duration is T1=3 minutes.
    # The maximal distance is D(3) = 3 + 1 = 4.
    # So, d(3) - d(0) = 4.
    T1 = 3
    dist_interval1 = T1 + 1
    d_3 = d_0 + dist_interval1

    # For the second interval [2, 7], the duration is T2=5 minutes.
    # The maximal distance is D(5) = 5 + 1 = 6.
    # So, d(7) - d(2) = 6.
    T2 = 5
    dist_interval2 = T2 + 1
    # This gives the equation: d(7) = d(2) + 6

    # Since the snail never turns back, the distance function d(t) is non-decreasing.
    # This means d(2) must be less than or equal to d(3).
    # d(2) <= d(3)
    # To maximize d(7), we must maximize d(2).
    # The maximum value for d(2) is d(3).
    # This corresponds to the snail not moving in the interval [2, 3].
    d_2_max = d_3

    # Now we calculate the maximal d(7).
    d_7_max = d_2_max + dist_interval2

    # Print out the logic and the final equation.
    print("Let d(t) be the distance at time t. We can set d(0) = 0.")
    print("Consider two overlapping time intervals: [0, 3] and [2, 7].")
    print("\nFor the interval [0, 3] (duration 3 minutes):")
    print("The maximum possible distance is 3 + 1 = 4 meters.")
    print(f"This means: d(3) - d(0) = {dist_interval1}")
    print(f"Since d(0) = 0, we have d(3) = {d_3}.")

    print("\nFor the interval [2, 7] (duration 5 minutes):")
    print("The maximum possible distance is 5 + 1 = 6 meters.")
    print(f"This means: d(7) - d(2) = {dist_interval2}")
    print(f"So, d(7) = d(2) + {dist_interval2}.")

    print("\nTo maximize d(7), we need to maximize d(2).")
    print("Since the snail never turns back, d(2) <= d(3).")
    print(f"The maximum value for d(2) is therefore d(3), which is {d_3}.")
    print(f"So, we set d(2) = {d_2_max}.")

    print("\nFinally, we can calculate the maximal distance at t=7:")
    print(f"d(7) = d(2) + {dist_interval2}")
    print(f"d(7) = {d_2_max} + {dist_interval2} = {d_7_max}")

    print(f"\nThe maximal distance is {d_7_max} meters.")

solve_snail_problem()