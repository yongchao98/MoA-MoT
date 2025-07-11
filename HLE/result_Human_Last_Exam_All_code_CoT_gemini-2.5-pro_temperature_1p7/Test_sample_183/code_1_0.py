import math
import random

# Caches to store results of computations to speed up the process,
# especially for primes that appear in the factorization of multiple numbers.
_prime_cache = {}
_inert_prime_cache = {}

def is_prime(n):
    """
    Miller-Rabin primality test.
    This is a probabilistic test, but for the given bases, it is deterministic
    for all integers up to a very high limit (well beyond 10^18).
    """
    if n in _prime_cache:
        return _prime_cache[n]
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    if n < 25:
        return True

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    
    # Bases for Miller-Rabin test sufficient for n up to ~3.3 * 10^24
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    for a in bases:
        if a >= n:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            _prime_cache[n] = False
            return False
            
    _prime_cache[n] = True
    return True

def pollard_rho(n):
    """
    Pollard's rho algorithm to find a non-trivial factor of a composite number n.
    It's particularly effective at finding small factors.
    """
    if n % 2 == 0:
        return 2
    
    # Use Floyd's cycle-finding algorithm with a random polynomial g(x) = (x^2 + c) mod n
    x = random.randint(1, n - 2)
    y = x
    c = random.randint(1, n - 1)
    g = 1
    f = lambda val: (pow(val, 2, n) + c) % n
    
    while g == 1:
        x = f(x)
        y = f(f(y))
        g = math.gcd(abs(x - y), n)

    if g == n:
        # Algorithm failed (found trivial factor), retry with different random values
        return pollard_rho(n)
    
    return g

def factorize(n):
    """
    Factorizes a large integer n into its prime factors and their exponents.
    Uses a stack-based iterative approach combining primality testing
    and Pollard's rho algorithm.
    """
    factors = {}
    num_stack = [n]
    
    while num_stack:
        current_num = num_stack.pop()
        
        if current_num == 1:
            continue

        if is_prime(current_num):
            factors[current_num] = factors.get(current_num, 0) + 1
            continue
            
        # If the number is composite, find a factor and add it and the
        # remaining part back to the stack to be factorized.
        factor = pollard_rho(current_num)
        num_stack.append(factor)
        num_stack.append(current_num // factor)
        
    return factors

def is_inert(p):
    """
    Checks if a prime p is inert in the number field Q(cuberoot(2)).
    A prime p is inert if p = 1 (mod 3) and 2 is not a cubic residue mod p.
    """
    if p in _inert_prime_cache:
        return _inert_prime_cache[p]
    
    # Primes p != 1 (mod 3) are not inert.
    if p % 3 != 1:
        _inert_prime_cache[p] = False
        return False
        
    # Use Euler's criterion for cubic residues: for a prime p = 1 (mod 3),
    # a is a cubic residue mod p iff a^((p-1)/3) = 1 (mod p).
    if pow(2, (p - 1) // 3, p) != 1:
        result = True
    else:
        result = False
        
    _inert_prime_cache[p] = result
    return result

# Main execution block
start_num = 10**18
end_num = 10**18 + 10000
representable_count = 0

for n in range(start_num, end_num + 1):
    prime_factors = factorize(n)
    is_representable = True
    for p, exponent in prime_factors.items():
        if is_inert(p):
            if exponent % 3 != 0:
                is_representable = False
                break
    if is_representable:
        representable_count += 1
        
print(f"Number of integers found: {representable_count}")