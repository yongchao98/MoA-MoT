import math

memo_p = {}

def partitions(n):
    """
    Calculates the number of partitions of an integer n using the recurrence
    from Euler's pentagonal number theorem. It uses memoization to store
    and retrieve previously computed values.
    """
    if n < 0:
        return 0
    if n == 0:
        return 1
    if n in memo_p:
        return memo_p[n]

    result = 0
    k = 1
    while True:
        pentagonal_1 = k * (3 * k - 1) // 2
        pentagonal_2 = k * (3 * k + 1) // 2
        
        # Stop if the next terms to subtract are for numbers larger than n
        if pentagonal_1 > n:
            break

        term1 = partitions(n - pentagonal_1)
        
        term2 = 0
        if pentagonal_2 <= n:
            term2 = partitions(n - pentagonal_2)

        if k % 2 != 0:  # If k is odd, add the terms
            result += term1 + term2
        else:  # If k is even, subtract the terms
            result -= term1 + term2
        
        k += 1
        
    memo_p[n] = result
    return result

def solve():
    """
    Solves the problem for n=9 by calculating the required values and printing the result.
    """
    n = 9
    
    # The number of conjugacy classes in Sigma_n is the partition number p(n).
    p_n = partitions(n)
    
    # The order of the symmetric group Sigma_n is n!
    fact_n = math.factorial(n)
    
    # The cardinality of Gamma(F) is (n!)^2 * p(n)
    cardinality = fact_n * fact_n * p_n
    
    print(f"For n = {n}:")
    print(f"The number of partitions, p({n}), which is the number of conjugacy classes of Sigma_{n}, is: {p_n}")
    print(f"The order of the symmetric group, {n}!, is: {fact_n}")
    print(f"The cardinality of [Gamma(F)](bullet,bullet) is given by the formula (n!)^2 * p(n).")
    print(f"Calculation: ({fact_n})^2 * {p_n} = {cardinality}")

solve()