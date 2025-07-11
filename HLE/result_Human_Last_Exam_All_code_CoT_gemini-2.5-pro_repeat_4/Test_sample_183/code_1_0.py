import random
import math

# This script may take several minutes to run due to the large number factorization.

def power(a, b, m):
    """Computes (a^b) % m using modular exponentiation."""
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = res * a
            if res >= m:
                res %= m
        a = a * a
        if a >= m:
            a %= m
        b //= 2
    return res

def is_prime_miller_rabin(n, k=5):
    """
    Miller-Rabin primality test.
    k is the number of rounds of testing.
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
        x = power(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = power(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def pollard_rho(n):
    """Pollard's rho algorithm to find a non-trivial factor of n."""
    if n % 2 == 0:
        return 2
    if is_prime_miller_rabin(n):
        return n
        
    x = random.randint(1, n - 2)
    y = x
    c = random.randint(1, n - 1)
    d = 1
    
    while d == 1:
        x = (power(x, 2, n) + c + n) % n
        y = (power(y, 2, n) + c + n) % n
        y = (power(y, 2, n) + c + n) % n
        d = math.gcd(abs(x - y), n)
        if d == n:
            # If d=n, the choice of x, y, c failed. Retry.
            return pollard_rho(n)
    return d

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    factors = {}
    
    # First, trial division for small primes
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            count = 0
            while n % p == 0:
                count += 1
                n //= p
            factors[p] = count

    # Use Pollard's rho for remaining factor
    def factorize_recursive(num, f):
        if num == 1:
            return
        if is_prime_miller_rabin(num):
            f[num] = f.get(num, 0) + 1
            return
        
        factor = pollard_rho(num)
        factorize_recursive(factor, f)
        factorize_recursive(num // factor, f)
        
    factorize_recursive(n, factors)
    return factors

memo_inert = {}
def is_inert(p):
    """Checks if a prime p is inert in Q(2^(1/3))."""
    if p in memo_inert:
        return memo_inert[p]
    
    # Primes 2 and 3 are not inert (they ramify).
    # Primes p = 2 (mod 3) are not inert (they split).
    if p % 3 != 1:
        memo_inert[p] = False
        return False
    
    # For p = 1 (mod 3), it's inert if 2 is not a cubic residue mod p.
    if power(2, (p - 1) // 3, p) != 1:
        memo_inert[p] = True
        return True
    
    memo_inert[p] = False
    return False

def solve():
    """
    Main function to solve the problem.
    """
    start_n = 10**18
    end_n = 10**18 + 10000
    
    count = 0
    
    for n in range(start_n, end_n + 1):
        factors = get_prime_factorization(n)
        
        is_representable = True
        for p, exponent in factors.items():
            if is_inert(p):
                if exponent % 3 != 0:
                    is_representable = False
                    break
        
        if is_representable:
            count += 1
            
    print(f"The number of such integers is: {count}")

solve()
<<<3334>>>