import math

def P(x):
    """
    Evaluates the polynomial P(X) at a given integer x.
    P(X) = (X^5-1)*(X^5-X)*(X^5-X^2)*(X^5-X^3)*(X^5-X^4)
    """
    if not isinstance(x, int):
        raise TypeError("Input must be an integer.")
    
    # Using python's arbitrary-precision integers
    term1 = x**5 - 1
    term2 = x**5 - x
    term3 = x**5 - x**2
    term4 = x**5 - x**3
    term5 = x**5 - x**4
    
    return term1 * term2 * term3 * term4 * term5

def is_prime(n):
    """
    Checks if a number is prime.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def get_prime_factorization_str(n):
    """
    Factorizes a number and returns a string representation like "2^10 * 3^2 * 5^1".
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    
    parts = []
    for p in sorted(factors.keys()):
        parts.append(f"{p}^{factors[p]}")
    return " * ".join(parts)

def solve():
    """
    Calculates the limit of the sequence g_n.
    """
    # We need to find the GCD of P(p) for all large primes p.
    # We can compute this numerically by taking the GCD over a range of primes.
    # Start with primes > 100.
    start_prime = 101
    
    # Find the first prime to initialize g
    current_p = start_prime
    while not is_prime(current_p):
        current_p += 1
    
    g = P(current_p)
    
    # Iterate through the next few primes and update g.
    # A range of primes up to 300 is sufficient to ensure the GCD stabilizes.
    current_p += 1
    end_p = 300 
    
    while current_p <= end_p:
        if is_prime(current_p):
            g = math.gcd(g, P(current_p))
        current_p += 1
        
    factor_str = get_prime_factorization_str(g)
    print(f"The limit is composed of the following prime powers:")
    print(f"{factor_str} = {g}")

solve()
<<<46080>>>