import math

def is_prime(n):
    """
    Miller-Rabin primality test. Sufficient for n up to ~3.3e23.
    """
    if n < 2:
        return False
    if n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
        
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
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
    Pollard's rho algorithm to find a non-trivial factor of n.
    """
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

        if d == n:
            # Factor not found, try a new c
            c += 1
            x, y = 2, 2
        elif d != 1:
            return d

def is_inert(p):
    """
    Checks if a prime p is inert in Z[cbrt(2)].
    Inert primes are p === 1 (mod 3) for which 2 is not a cubic residue mod p.
    """
    if p <= 3: # 2 and 3 are not inert
        return False
    if p % 3 != 1:
        return False
    if pow(2, (p - 1) // 3, p) != 1:
        return True
    return False

def get_prime_factors(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    d = [n]
    while d:
        c = d.pop()
        if c == 1:
            continue
        if is_prime(c):
            factors[c] = factors.get(c, 0) + 1
            continue
        f = pollard_rho(c)
        d.append(f)
        d.append(c // f)
    return factors

def is_representable(n, primes_list):
    """
    Checks if n can be expressed in the form x^3 + 2y^3 + 4z^3 - 6xyz.
    """
    temp_n = n

    # Trial division for small primes
    for p in primes_list:
        if p * p > temp_n:
            break
        if temp_n % p == 0:
            count = 0
            while temp_n % p == 0:
                count += 1
                temp_n //= p
            if is_inert(p) and count % 3 != 0:
                return False

    # Handle the remaining large factor m
    m = temp_n
    if m == 1:
        return True

    # Get prime factors of the remaining number m
    factors_m = get_prime_factors(m)
    for p, exponent in factors_m.items():
        if is_inert(p) and exponent % 3 != 0:
            return False

    return True

def solve():
    """
    Main function to solve the problem.
    """
    LOWER_BOUND = 10**18
    UPPER_BOUND = 10**18 + 10000
    SIEVE_LIMIT = 1000000  # Primes up to 10^6 for trial division

    # Generate primes using a sieve
    primes_list = []
    is_prime_sieve = [True] * (SIEVE_LIMIT + 1)
    is_prime_sieve[0] = is_prime_sieve[1] = False
    for p in range(2, SIEVE_LIMIT + 1):
        if is_prime_sieve[p]:
            primes_list.append(p)
            for i in range(p * p, SIEVE_LIMIT + 1, p):
                is_prime_sieve[i] = False

    count = 0
    for n in range(LOWER_BOUND, UPPER_BOUND + 1):
        if is_representable(n, primes_list):
            count += 1
    
    # As per instructions, "output each number in the final equation"
    # which here is the total count.
    print(count)

if __name__ == '__main__':
    solve()
<<<1113>>>