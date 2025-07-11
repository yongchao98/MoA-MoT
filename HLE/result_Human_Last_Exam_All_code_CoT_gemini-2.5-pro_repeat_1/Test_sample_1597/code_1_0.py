def count_divisors(n):
    """
    Computes the number of divisors of a positive integer n, denoted as tau(n).
    """
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    num_divs = 1
    temp_n = n
    
    # Efficiently find prime factorization and calculate tau(n).
    # tau(p1^a1 * p2^a2 * ...) = (a1+1)*(a2+1)*...
    
    # Handle the factor 2
    power = 0
    while temp_n % 2 == 0:
        power += 1
        temp_n //= 2
    if power > 0:
        num_divs *= (power + 1)
        
    # Handle odd factors
    d = 3
    while d * d <= temp_n:
        if temp_n % d == 0:
            power = 0
            while temp_n % d == 0:
                power += 1
                temp_n //= d
            num_divs *= (power + 1)
        d += 2
        
    # If temp_n is still greater than 1, it must be a prime factor itself.
    if temp_n > 1:
        num_divs *= 2
        
    return num_divs

def solve():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    """
    N = 1000
    
    # Calculate Sum(tau(k)) for k from 1 to N
    sum_of_tau = 0
    for k in range(1, N + 1):
        sum_of_tau += count_divisors(k)
    
    # The total number of allowed pairs is 2 * |S1| - |S1 intersect S2|
    # which is 2 * sum_of_tau - N
    result = 2 * sum_of_tau - N
    
    # Output the numbers in the final equation
    print(f"The total number of divisors sum, Sum(tau(k) for k=1..1000), is: {sum_of_tau}")
    print(f"The number of allowed pairs is 2 * {sum_of_tau} - {N} = {result}")

solve()