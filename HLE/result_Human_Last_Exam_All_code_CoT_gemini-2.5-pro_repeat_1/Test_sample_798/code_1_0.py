import math

def solve_partition_problem():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    # Step 1: Define the range and calculate the total number of integers (N).
    start = 10001
    end = 42149572
    num_integers = end - start + 1

    # Step 2: Determine the maximum possible size (m_max) of a single set.
    # The maximum possible distance (D) in the range is end - start.
    max_distance = end - start

    # A set of size 'm' has m*(m-1)/2 distinct distances. These distances must be
    # less than or equal to the maximum possible distance in the range.
    # So, we find the largest integer 'm' that satisfies: m * (m - 1) / 2 <= max_distance
    # This is equivalent to solving the quadratic equation m^2 - m - 2*max_distance = 0
    # and taking the floor of the positive root.
    # m = (1 + sqrt(1 + 8 * max_distance)) / 2
    max_set_size = math.floor((1 + math.sqrt(1 + 8 * max_distance)) / 2)

    # Step 3: Calculate the minimum number of sets (k).
    # This is the total number of integers divided by the maximum possible size of one set,
    # rounded up to the nearest whole number (ceiling).
    min_sets = math.ceil(num_integers / max_set_size)

    # Print the logic and the final equation.
    print("Problem: Find the minimum number of distinct-distance-sets to partition the integers from 10001 to 42149572.")
    print("\nStep 1: Calculate the total number of integers (N).")
    print(f"N = {end} - {start} + 1 = {num_integers}")

    print("\nStep 2: Find the maximum possible size (m_max) for a single set.")
    print(f"The maximum distance (D) is {end} - {start} = {max_distance}.")
    print("The number of distances in a set of size m is m*(m-1)/2. This must be <= D.")
    print(f"Solving m*(m-1)/2 <= {max_distance} for the largest integer m gives:")
    print(f"m_max = {max_set_size}")

    print("\nStep 3: Calculate the minimum number of sets (k).")
    print("This is the lower bound found by dividing N by m_max and taking the ceiling.")
    print("Final Equation:")
    print(f"k = ceil(N / m_max)")
    print(f"k = ceil({num_integers} / {max_set_size})")
    print(f"k = {int(min_sets)}")

solve_partition_problem()
<<<4591>>>