def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of the complement of X in SO(4)
    for degrees up to 100.
    """
    total_rank = 0
    # The Poincare series for the cohomology ring is (t^3 + t^6) / (1 - t^4)^2
    # We expand this and sum the coefficients for degrees k from 0 to 100.
    # The expansion of 1/(1-t^4)^2 is sum_{n=0 to inf} (n+1)t^{4n}.
    
    # Contribution from t^3 * sum_{n=0 to inf} (n+1)t^{4n} = sum_{n=0 to inf} (n+1)t^{4n+3}
    # We need to find the sum of ranks for degrees k <= 100.
    # k = 4n+3 <= 100  => 4n <= 97 => n <= 24.25. So, n goes from 0 to 24.
    # The rank for a given n is n+1.
    sum1 = 0
    for n in range(25):
        k = 4 * n + 3
        rank = n + 1
        # print(f"Degree k={k}, rank={rank}") # for debugging
        sum1 += rank
        
    # Contribution from t^6 * sum_{n=0 to inf} (n+1)t^{4n} = sum_{n=0 to inf} (n+1)t^{4n+6}
    # We need to find the sum of ranks for degrees k <= 100.
    # k = 4n+6 <= 100 => 4n <= 94 => n <= 23.5. So, n goes from 0 to 23.
    # The rank for a given n is n+1.
    sum2 = 0
    for n in range(24):
        k = 4 * n + 6
        rank = n + 1
        # print(f"Degree k={k}, rank={rank}") # for debugging
        sum2 += rank
        
    total_rank = sum1 + sum2
    print(f"The sum of ranks from the t^3 term is: {sum1}")
    print(f"The sum of ranks from the t^6 term is: {sum2}")
    print(f"The total rank of A in degrees * <= 100 is: {total_rank}")

solve()