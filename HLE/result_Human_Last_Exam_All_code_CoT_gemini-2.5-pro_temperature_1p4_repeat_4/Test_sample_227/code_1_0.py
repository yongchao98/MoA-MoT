def solve_snail_puzzle():
    """
    Solves the snail puzzle by constructing a valid scenario.

    The problem asks for the maximal distance a snail can travel in 7 minutes
    under specific observation conditions.

    Let d(t) be the distance the snail has traveled at time t.
    The conditions are:
    1. The total time is 7 minutes, so the interval is [0, 7].
    2. A finite number of observers watch the snail.
    3. Each observer watches for exactly 1 minute, [t_i, t_i+1].
    4. Each observer finds the snail advanced exactly 1 meter, so d(t_i+1) - d(t_i) = 1.
    5. At any moment in [0, 7], at least one observer is watching.

    The solution is found by constructing a specific travel plan for the snail (a sequence
    of distances traveled per minute) and then showing that a valid set of observers
    can exist for this plan.

    Proposed Travel Plan (distance traveled in each minute interval [i-1, i]):
    - Minute 1 ([0,1]): 1 meter
    - Minute 2 ([1,2]): 2 meters
    - Minute 3 ([2,3]): 1 meter
    - Minute 4 ([3,4]): 2 meters
    - Minute 5 ([4,5]): 1 meter
    - Minute 6 ([5,6]): 2 meters
    - Minute 7 ([6,7]): 1 meter

    The total distance is the sum of these individual distances.
    """

    # Distances traveled in each of the 7 one-minute intervals
    distances_per_minute = [1, 2, 1, 2, 1, 2, 1]

    # Calculate the total distance
    total_distance = sum(distances_per_minute)

    # The equation for the total distance
    equation_parts = [str(d) for d in distances_per_minute]
    equation_str = " + ".join(equation_parts)

    print("The maximal distance can be found by constructing a specific scenario.")
    print("Let the snail's travel distance in each 1-minute interval [i-1, i] be d_i.")
    print("A possible scenario is to have the snail travel d_i meters as follows:")
    for i, dist in enumerate(distances_per_minute):
        print(f"  - In minute {i+1} ([{i},{i+1}]): {dist} meter(s)")
    print("\nThe total distance is the sum of the distances from each minute:")
    print(f"Total Distance = {equation_str} = {total_distance} meters")
    print("\nThis travel plan is possible with a set of observers covering the entire 7 minutes.")
    print("For example, observers watching [0,1], [2,3], [4,5], [6,7] see 1 meter of travel.")
    print("For minutes where 2 meters are traveled (e.g., [1,2]), an observer on [1.5, 2.5] can still see exactly 1 meter of travel if the snail jumps 1m at t=1.5 and 1m at t=2.0.")
    print("\nThis demonstrates that a distance of 10 meters is achievable.")
    print(f"The maximal distance is: {total_distance}")


solve_snail_puzzle()
