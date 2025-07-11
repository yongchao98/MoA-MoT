import math

def is_prime(n, k=5):
    """Miller-Rabin primality test."""
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

    # Pre-selected bases for numbers up to 10^18
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n < 341550071728321: # This covers numbers up to ~3.4e14
        bases = [2, 3, 5, 7, 11, 13, 17]
    if n < 4759123141: # This covers numbers up to ~4.7e9
       bases = [2, 7, 61]


    for a in bases:
        if a >= n: break
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
    """Pollard's rho algorithm to find a non-trivial factor."""
    if n % 2 == 0:
        return 2
    if is_prime(n):
        return n
    
    x = 2
    y = 2
    c = 1
    d = 1
    while True:
        x = (pow(x, 2, n) + c) % n
        y = (pow(y, 2, n) + c) % n
        y = (pow(y, 2, n) + c) % n
        d = math.gcd(abs(x - y), n)
        if d == 1:
            continue
        if d != n:
            return d
        # d == n, failure, try another c
        c += 1
        x, y = 2, 2


memo_factors = {}
def has_no_2mod3_factors(m):
    """Check if m has any prime factors q = 2 (mod 3)."""
    if m in memo_factors:
        return memo_factors[m]
    
    if m == 1:
        return True
    if m % 3 == 2:
        memo_factors[m] = False
        return False

    s = math.isqrt(m)
    if s * s == m:
        res = has_no_2mod3_factors(s)
        memo_factors[m] = res
        return res

    if is_prime(m):
        # m % 3 is 1, so it's a valid prime factor
        memo_factors[m] = True
        return True

    factor = pollard_rho(m)
    res = has_no_2mod3_factors(factor) and has_no_2mod3_factors(m // factor)
    memo_factors[m] = res
    return res

def solve():
    """Main function to solve the problem."""
    
    # Generate primes up to 10^6
    limit = 1000000
    primes = []
    is_prime_sieve = [True] * (limit + 1)
    is_prime_sieve[0] = is_prime_sieve[1] = False
    for p in range(2, limit + 1):
        if is_prime_sieve[p]:
            primes.append(p)
            for i in range(p * p, limit + 1, p):
                is_prime_sieve[i] = False

    count = 0
    start_n = 10**18
    end_n = 10**18 + 10000

    for k in range(end_n - start_n + 1):
        n = start_n + k
        
        m = n
        is_rep = True

        for p in primes:
            if p * p > m:
                break
            if m % p == 0:
                exp = 0
                while m % p == 0:
                    exp += 1
                    m //= p
                if p % 3 == 2 and exp % 3 != 0:
                    is_rep = False
                    break
        
        if not is_rep:
            continue
            
        # After dividing by small primes, check the remainder m
        # As argued, m cannot have prime factors q > 10^6 with q=2(mod 3)
        if m > 1:
            if not has_no_2mod3_factors(m):
                is_rep = False

        if is_rep:
            count += 1
            
    print(f"The number of such integers is: {count}")

solve()
<<<3334>>>