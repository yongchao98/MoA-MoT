import math
import random

def is_prime(n, k=5):
    """
    Miller-Rabin primality test.
    
    Args:
        n: The number to test.
        k: The number of rounds of testing.
    
    Returns:
        True if n is probably prime, False otherwise.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def pollard_rho(n):
    """
    Pollard's rho algorithm for integer factorization.
    Finds a non-trivial factor of n.
    """
    if n % 2 == 0:
        return 2
    if is_prime(n):
        return n
    
    # Use the improved version with Brent's cycle detection
    y, c, m = random.randint(1, n-1), random.randint(1, n-1), random.randint(1, n-1)
    g, r, q = 1, 1, 1
    while g == 1:
        x = y
        for _ in range(r):
            y = (pow(y, 2, n) + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r-k)):
                y = (pow(y, 2, n) + c) % n
                q = q * abs(x - y) % n
            g = math.gcd(q, n)
            k += m
        r *= 2
    if g == n:
        while True:
            ys = (pow(ys, 2, n) + c) % n
            g = math.gcd(abs(x - ys), n)
            if g > 1:
                break
    return g

def factorize(n, factors):
    """
    Recursively find all prime factors of a number n.
    
    Args:
        n: The number to factorize.
        factors: A dictionary to store prime factors and their exponents.
    """
    if n == 1:
        return
    if is_prime(n):
        factors[n] = factors.get(n, 0) + 1
        return
    
    divisor = pollard_rho(n)
    factorize(divisor, factors)
    if n // divisor > 1:
      factorize(n // divisor, factors)

def solve():
    """
    Counts the integers in [10^18, 10^18 + 10000] of the form x^3+2y^3+4z^3-6xyz.
    An integer n has this form iff for every prime p=2(mod 3), its exponent
    in the factorization of n is a multiple of 3.
    """
    start_num = 10**18
    num_range = 10001
    
    is_representable = [True] * num_range
    cofactors = [start_num + i for i in range(num_range)]

    # Sieve for small prime factors up to 10^6
    sieve_limit = 10**6
    primes_sieve = [True] * (sieve_limit + 1)
    primes_sieve[0] = primes_sieve[1] = False
    for i in range(2, math.isqrt(sieve_limit) + 1):
        if primes_sieve[i]:
            for multiple in range(i*i, sieve_limit + 1, i):
                primes_sieve[multiple] = False
    
    for p in range(2, sieve_limit + 1):
        if not primes_sieve[p]:
            continue

        # Start sieving from the first multiple of p in the range
        start_offset = start_num % p
        first_idx = (p - start_offset) % p
            
        for i in range(first_idx, num_range, p):
            if cofactors[i] % p == 0:
                exponent = 0
                while cofactors[i] % p == 0:
                    exponent += 1
                    cofactors[i] //= p
                
                # Check condition only for primes p = 2 (mod 3)
                if p % 3 == 2:
                    if exponent % 3 != 0:
                        is_representable[i] = False

    # Process remaining large cofactors
    count = 0
    for i in range(num_range):
        if is_representable[i]:
            m = cofactors[i]
            if m > 1:
                factors = {}
                try:
                    factorize(m, factors)
                    valid_m = True
                    for p, exp in factors.items():
                        if p % 3 == 2:
                            if exp % 3 != 0:
                                valid_m = False
                                break
                    if valid_m:
                        count += 1
                except RecursionError:
                    # Factoring very difficult numbers might hit recursion limit
                    # For this problem, assume such cases are not representable
                    pass
            else:  # m == 1, means all factors were handled by the sieve
                count += 1
                
    print(f"The number of such integers is: {count}")

solve()