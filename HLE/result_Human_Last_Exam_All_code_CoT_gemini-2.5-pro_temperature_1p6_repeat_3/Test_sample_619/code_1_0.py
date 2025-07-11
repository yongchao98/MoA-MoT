import math

def solve():
    """
    Solves the problem by calculating the minimum value of the sum of sizes of the sets.
    The solution is based on the interpretation that the sets are intervals on the real line
    and their size is their measure.
    """
    try:
        n_str = input("Enter the value of n (an integer >= 1): ")
        n = int(n_str)
        if n < 1:
            print("Please enter an integer greater than or equal to 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The problem is to find the minimum of sum(|S_i|) given |S_i △ S_j| = |i-j|.
    # This can be modeled by setting S_i = [min(i,k), max(i,k)] for a pivot k.
    # The sum is then sum(|i-k|) for i=1 to n.
    # This sum is minimized when k is the median of {1, ..., n}.
    # The minimum value of this sum is equal to floor(n^2 / 4).
    
    min_sum = math.floor(n**2 / 4)

    if n % 2 == 0:
        # n = 2m
        m = n // 2
        # k can be m or m+1. Let's use k=m.
        # The sum is sum_{i=1 to 2m} |i-m|.
        # sum from i=1 to m-1 is (m-1) + (m-2) + ... + 1 = m(m-1)/2
        # sum from i=m+1 to 2m is 1 + 2 + ... + m = m(m+1)/2
        # Total sum is m(m-1)/2 + 0 + m(m+1)/2 = m^2
        print(f"For n = {n}, which is even, we can write n = 2 * {m}.")
        print(f"The minimum sum is given by the formula m^2.")
        print(f"Calculation: {m}^2 = {m*m}")
    else:
        # n = 2m+1
        m = (n - 1) // 2
        # The median k is m+1.
        # The sum is sum_{i=1 to 2m+1} |i-(m+1)|.
        # This is 2 * (1 + 2 + ... + m) = 2 * m(m+1)/2 = m(m+1).
        print(f"For n = {n}, which is odd, we can write n = 2 * {m} + 1.")
        print(f"The minimum sum is given by the formula m*(m+1).")
        print(f"Calculation: {m} * ({m}+1) = {m} * {m+1} = {m*(m+1)}")
    
    print("\nThis result matches the formula floor(n^2 / 4):")
    print(f"floor({n}^2 / 4) = floor({n*n} / 4) = floor({n*n/4}) = {min_sum}")
    print(f"\nThe minimum value of the sum is: {min_sum}")

solve()