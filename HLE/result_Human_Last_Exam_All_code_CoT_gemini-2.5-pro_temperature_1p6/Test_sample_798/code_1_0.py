import math

def solve_partition_problem():
    """
    Calculates the minimum number of distinct distance sets needed to partition
    the integers from 10001 to 42149572.
    """
    
    # The range of integers to be partitioned
    N1 = 10001
    N2 = 42149572

    # Step 1: Calculate the total number of integers in the set.
    num_integers = N2 - N1 + 1

    # Step 2: Determine the maximum possible size of a single distinct distance set.
    # A "distinct distance set" S of size n is a set of integers {s_1, ..., s_n}
    # where all distances |s_i - s_j| for i != j are unique.
    # For a sorted set s_1 < ... < s_n, there are n*(n-1)/2 distinct distances.
    # The largest distance, s_n - s_1, must be at least n*(n-1)/2.
    # The elements are from [N1, N2], so the maximum possible span for any
    # subset is N2 - N1.
    span = N2 - N1

    # This gives us the inequality for the maximum possible size 'n' of a set:
    # n*(n-1)/2 <= span
    # n^2 - n - 2*span <= 0
    
    # To find the maximum integer 'n' (n_max) satisfying this, we solve the 
    # quadratic equation x^2 - x - 2*span = 0. The integer part of the positive 
    # root gives us n_max.
    c = 2 * span
    # Positive root of x^2 - x - c = 0 is (1 + sqrt(1 + 4c)) / 2
    root = (1 + math.sqrt(1 + 4 * c)) / 2
    n_max = int(root)

    # Verification to ensure we have the largest integer n satisfying the condition
    if (n_max + 1) * n_max <= c:
        n_max += 1

    # Step 3: Calculate the minimum number of sets.
    # To partition a set of 'num_integers' into subsets of maximum size 'n_max',
    # we need at least ceil(num_integers / n_max) subsets.
    min_sets = math.ceil(num_integers / n_max)

    print("--- Problem Analysis ---")
    print(f"The integer range to partition is [{N1}, {N2}].")
    print(f"The total number of integers is {N2} - {N1} + 1 = {num_integers}.")
    print(f"The span of the interval is {N2} - {N1} = {span}.")
    print("\n--- Calculating Maximum Set Size ---")
    print(f"The size 'n' of any distinct distance set must satisfy: n*(n-1)/2 <= {span}")
    print(f"This is equivalent to: n*(n-1) <= {c}")
    print(f"The maximum integer 'n' that satisfies this inequality is: {n_max}")

    print("\n--- Calculating Minimum Number of Sets ---")
    print("The minimum number of sets is the ceiling of the total number of integers divided by this maximum size.")
    print(f"Final calculation: ceil({num_integers} / {n_max}) = {min_sets}")
    
solve_partition_problem()