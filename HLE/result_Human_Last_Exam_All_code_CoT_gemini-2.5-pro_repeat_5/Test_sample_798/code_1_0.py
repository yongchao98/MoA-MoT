import math

def solve_partition_problem():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Calculate N, the total number of integers to partition.
    num_integers_N = end - start + 1

    # Calculate L, the length of the interval, which is the maximum possible range for any subset.
    interval_length_L = end - start

    # Find the maximum possible size (k_max) of a DDS.
    # We need to find the largest integer k such that: k * (k - 1) / 2 <= L
    # This is equivalent to solving for the positive root of the quadratic equation:
    # k^2 - k - 2*L = 0
    # The root is (1 + sqrt(1 + 8*L)) / 2. We take the floor to get the max integer k.
    discriminant = 1 + 8 * interval_length_L
    k_max_float = (1 + math.sqrt(discriminant)) / 2
    k_max = math.floor(k_max_float)

    # Calculate the minimum number of sets (m) required for the partition.
    # m = ceil(N / k_max)
    min_sets_m = math.ceil(num_integers_N / k_max)

    print("--- Problem Analysis ---")
    print(f"Total integers to partition (N): {num_integers_N}")
    print(f"Maximum possible range for a set (L): {interval_length_L}")
    print(f"Maximum possible size of a single distinct-distance-set (k_max): {k_max}")
    print("\n--- Final Calculation ---")
    print(f"The minimum number of sets is calculated as ceil(N / k_max).")
    print(f"Equation: ceil({num_integers_N} / {k_max}) = {min_sets_m}")
    print(f"\nThe minimum number of distinct-distance-sets needed is {min_sets_m}.")

solve_partition_problem()