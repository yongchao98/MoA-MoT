import math

def solve_partition_problem():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1: Calculate the total number of integers (N) and the maximum span (L).
    N = end - start + 1
    L = end - start

    print(f"The range of integers is from {start} to {end}.")
    print(f"Total number of integers to partition (N) = {end} - {start} + 1 = {N}")
    print(f"Maximum possible span for a set (L) = {end} - {start} = {L}")
    print("-" * 20)

    # Step 2: Find the maximum possible size (k_max) of a distinct-distance-set (DDS).
    # The number of distances in a set of size k is k*(k-1)/2.
    # This must be less than or equal to the span L.
    # k_max * (k_max - 1) / 2 <= L
    # k_max^2 - k_max - 2*L <= 0
    # We solve for k_max using the quadratic formula for the roots of x^2 - x - 2L = 0,
    # which are (1 +/- sqrt(1 + 8*L)) / 2. We take the positive root.
    
    discriminant = 1 + 8 * L
    k_max_float = (1 + math.sqrt(discriminant)) / 2
    k_max = math.floor(k_max_float)

    print("To find the maximum size of a DDS (k_max), we use the inequality:")
    print("k_max * (k_max - 1) / 2 <= L")
    print(f"k_max * (k_max - 1) <= 2 * {L}")
    print(f"k_max * (k_max - 1) <= {2 * L}")
    print(f"Solving for k_max gives a value of approximately {k_max_float:.4f}.")
    print(f"The maximum integer size for a set is k_max = floor({k_max_float:.4f}) = {k_max}")
    print("-" * 20)

    # Step 3: Calculate the minimum number of sets (m).
    # To partition N integers into sets of maximum size k_max, we need at least N / k_max sets.
    # Since the number of sets must be an integer, we take the ceiling.
    
    min_sets_float = N / k_max
    min_sets = math.ceil(min_sets_float)

    print("The minimum number of sets (m) is found by dividing the total number of integers (N)")
    print("by the maximum possible size of a single set (k_max).")
    print("m >= N / k_max")
    print(f"m >= {N} / {k_max}")
    print(f"m >= {min_sets_float:.4f}")
    print(f"Since the number of sets must be an integer, we take the ceiling.")
    print(f"Minimum number of sets = ceil({min_sets_float:.4f}) = {min_sets}")
    print("-" * 20)
    
    print(f"Final Answer: The minimum number of distinct-distance-sets required is {min_sets}.")

solve_partition_problem()