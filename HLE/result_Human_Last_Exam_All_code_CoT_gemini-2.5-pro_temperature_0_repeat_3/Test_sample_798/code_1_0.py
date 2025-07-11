import math

def solve():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1 & 2: Calculate the total number of integers, N.
    # The problem is equivalent to partitioning the set {1, 2, ..., N}.
    N = end - start + 1
    print(f"The range of integers is from {start} to {end}.")
    print(f"The total number of integers to partition is N = {end} - {start} + 1 = {N}.")
    print("-" * 20)

    # Step 3 & 4: Check if N is of the form m * (m - 1).
    # This requires solving the quadratic equation m^2 - m - N = 0 for an integer m.
    # The solution for m is (1 + sqrt(1 + 4*N)) / 2.
    discriminant = 1 + 4 * N
    sqrt_discriminant = math.isqrt(discriminant)

    # Step 5: If N fits the form, calculate m.
    if sqrt_discriminant * sqrt_discriminant == discriminant:
        # Check if (1 + sqrt_discriminant) is even to ensure m is an integer.
        if (1 + sqrt_discriminant) % 2 == 0:
            m = (1 + sqrt_discriminant) // 2
            print(f"The number N can be expressed in the form m * (m - 1).")
            print(f"We found that for m = {m}, m * (m - 1) = {m} * {m-1} = {N}.")
            print("-" * 20)

            # Step 6: Apply the known result/conjecture for k(N) when N = m*(m-1).
            # The minimum number of Sidon sets required to partition {1, ..., m*(m-1)}
            # is k = m - 1.
            k = m - 1
            
            print(f"The minimum number of distinct-distance-sets is given by the formula k = m - 1.")
            print(f"The final equation is:")
            print(f"{m} - 1 = {k}")
            print("-" * 20)
            print(f"The final answer is {k}.")

        else:
            print("Error: N fits the pattern, but m is not an integer.")
    else:
        print("Error: N is not of the form m*(m-1).")

solve()