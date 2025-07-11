import math
import random

def is_prime(n, k=5):
    """
    Miller-Rabin primality test.
    
    Args:
        n: The number to test.
        k: The number of rounds of testing.
        
    Returns:
        True if n is likely prime, False otherwise.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
        
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def pollard_rho(n):
    """
    Pollard's rho algorithm to find a non-trivial factor.
    Assumes n is composite.
    
    Args:
        n: The number to factor.
        
    Returns:
        A non-trivial factor of n.
    """
    if n % 2 == 0:
        return 2
    
    # Use a non-deterministic version for robustness
    while True:
        x = random.randint(1, n - 2)
        y = x
        c = random.randint(1, n - 1)
        d = 1
        
        def f(val):
            return (pow(val, 2, n) + c) % n

        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
        
        if d != n:
            return d

def factorize_range(start, limit, primes_list):
    """
    Factorizes numbers in a given range using a sieve followed by
    factorization of large remainders.
    
    Args:
        start: The starting number of the range.
        limit: The number of integers in the range.
        primes_list: A list of primes to use for sieving.
        
    Returns:
        A list of dictionaries, where each dictionary holds the prime
        factorization of the corresponding number in the range.
    """
    numbers = [start + i for i in range(limit)]
    factors = [{} for _ in range(limit)]

    for p in primes_list:
        # Find the first multiple of p in the range
        start_offset = (start + p - 1) // p
        start_idx = start_offset * p - start
        
        for i in range(start_idx, limit, p):
            e = 0
            while numbers[i] % p == 0:
                numbers[i] //= p
                e += 1
            if e > 0:
                factors[i][p] = e

    for i in range(limit):
        m = numbers[i]
        if m > 1:
            # Process the remaining number m, which has only large prime factors
            if is_prime(m):
                factors[i][m] = factors[i].get(m, 0) + 1
            else:
                s = math.isqrt(m)
                if s * s == m:
                    # m is a perfect square of a large prime
                    factors[i][s] = factors[i].get(s, 0) + 2
                else:
                    # m is a product of two large distinct primes
                    p = pollard_rho(m)
                    q = m // p
                    factors[i][p] = factors[i].get(p, 0) + 1
                    factors[i][q] = factors[i].get(q, 0) + 1
    
    return factors

def is_inert_prime(p):
    """
    Checks if a prime p is inert in the number field Q(2^(1/3)).
    Condition: p = 1 (mod 3) and 2 is not a cubic residue mod p.
    """
    if p <= 3:  # Ramified or special primes
        return False
    if p % 3 != 1:
        return False
    
    # Check if 2 is a cubic residue mod p
    if pow(2, (p - 1) // 3, p) != 1:
        return True
        
    return False

def solve():
    """
    Main function to solve the problem.
    """
    start = 10**18
    limit = 10001

    # Generate primes for sieving. A bound of 2*10^6 is sufficient to ensure
    # that the remainder has at most two prime factors.
    sieve_bound = 2 * 10**6 
    primes = []
    is_p = [True] * (sieve_bound + 1)
    is_p[0] = is_p[1] = False
    for p in range(2, sieve_bound + 1):
        if is_p[p]:
            primes.append(p)
            for i in range(p * p, sieve_bound + 1, p):
                is_p[i] = False
    
    # Get prime factorizations for all numbers in the range
    all_factors = factorize_range(start, limit, primes)
    
    count = 0
    
    for i in range(limit):
        is_a_norm = True
        factors = all_factors[i]
        for p, e in factors.items():
            if is_inert_prime(p):
                if e % 3 != 0:
                    is_a_norm = False
                    break
        if is_a_norm:
            count += 1
            
    print(f"Searching in range: [{start}, {start + limit - 1}]")
    print(f"Total numbers checked: {limit}")
    print(f"Determined how many integers can be expressed in the form n = x^3 + 2y^3 + 4z^3 - 6xyz.")
    print(f"Result: {count}")
    
    return count

# Run the solver
final_count = solve()
# The final answer format as requested by the prompt.
# <<<final_count>>>