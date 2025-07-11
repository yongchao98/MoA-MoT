import math

def calculate_min_sets():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start_num = 10001
    end_num = 42149572

    # The total number of integers to partition.
    total_integers = end_num - start_num + 1

    # The maximum possible range (span) for any subset of these integers.
    max_span = end_num - start_num

    # For a distinct distance set of size 'n', the span of the set must be at least n*(n-1)/2.
    # Therefore, for any set in our partition, its size 'n' must satisfy:
    # n * (n - 1) / 2 <= max_span
    # We need to find the maximum integer 'n' (n_max) that satisfies this inequality.
    # The inequality can be written as: n^2 - n - 2 * max_span <= 0.
    # We find the positive root of n^2 - n - 2*max_span = 0 to get the upper bound for n.
    discriminant = 1 + 8 * max_span
    n_upper_bound = (1 + math.sqrt(discriminant)) / 2
    n_max = math.floor(n_upper_bound)

    # With the total number of integers (total_integers) and the maximum size of any set (n_max),
    # the minimum number of sets required is the ceiling of their division.
    min_num_sets_float = total_integers / n_max
    min_num_sets = math.ceil(min_num_sets_float)
    
    print("Step 1: Define the range and counts.")
    print(f"Integers are partitioned from {start_num} to {end_num}.")
    print(f"Total number of integers (N) = {end_num} - {start_num} + 1 = {total_integers}")
    print(f"Maximum span for any set (L) = {end_num} - {start_num} = {max_span}")
    print("-" * 50)

    print("Step 2: Find the maximum possible size for a single set (n_max).")
    print("The constraint is: n * (n - 1) / 2 <= L")
    print(f"n * (n - 1) <= 2 * {max_span}")
    print(f"n_max = floor((1 + sqrt(1 + 8 * {max_span})) / 2) = {n_max}")
    print("-" * 50)

    print("Step 3: Calculate the minimum number of sets.")
    print("The minimum number of sets is ceil(N / n_max).")
    print(f"min_sets = ceil({total_integers} / {n_max})")
    print(f"min_sets = ceil({min_num_sets_float})")
    print(f"min_sets = {min_num_sets}")
    print("-" * 50)

    print(f"The minimum number of distinct-distance-sets required is {min_num_sets}.")


calculate_min_sets()