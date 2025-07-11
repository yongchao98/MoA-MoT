import math

def solve():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1: Calculate the total number of integers to partition.
    N = end - start + 1
    print(f"The range of integers is from {start} to {end}.")
    print(f"The total number of integers to partition is N = {end} - {start} + 1 = {N}.")
    print("-" * 20)

    # Step 2: Check if N is of the form m * (m-1) by solving m^2 - m - N = 0.
    # The positive solution for m is (1 + sqrt(1 + 4*N)) / 2.
    delta = 1 + 4 * N
    sqrt_delta = math.isqrt(delta)

    # Check if delta is a perfect square.
    if sqrt_delta * sqrt_delta != delta:
        print("N is not of the form m*(m-1) for an integer m.")
        # If not, the problem might be more complex, but we proceed with approximation.
        # k >= sqrt(N)
        k_min = math.ceil(math.sqrt(N))
        print(f"A lower bound for the number of sets is ceil(sqrt(N)) = {k_min}.")
        return

    m = (1 + sqrt_delta) // 2
    
    print("The number of integers N can be expressed in the form m * (m-1).")
    print(f"We solve the equation: m * (m - 1) = {N}")
    print(f"This is equivalent to the quadratic equation: m^2 - m - {N} = 0")
    print(f"The positive integer solution for m is {m}.")
    print(f"Verification: {m} * ({m} - 1) = {m * (m-1)}")
    print("-" * 20)

    # Step 3: Based on the pattern for N = m*(m-1), the minimum number of sets is m-1.
    result = m - 1
    print("For N = m*(m-1), the minimum number of distinct-distance-sets (or Sidon sets) required for a partition is m-1.")
    print(f"The final answer is m - 1 = {m} - 1 = {result}.")

solve()
