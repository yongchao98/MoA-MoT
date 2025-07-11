import math

# --- Number Theory Utility Functions ---

_known_primes_for_mr = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

def is_prime(n):
    """
    Miller-Rabin primality test for 64-bit integers.
    """
    if n < 2:
        return False
    if n in _known_primes_for_mr:
        return True
    if any(n % p == 0 for p in _known_primes_for_mr):
        return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for p in _known_primes_for_mr:
        if p >= n:
            break
        x = pow(p, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

_memo_is_inert = {}
def is_inert(p):
    """
    Checks if a prime p is inert in the number field Q(2^(1/3)).
    Memoization is used to cache results.
    """
    if p in _memo_is_inert:
        return _memo_is_inert[p]
    
    if not is_prime(p):
        _memo_is_inert[p] = False
        return False
    
    # Primes p such that t^3-2 is irreducible mod p
    # are those where p=1(mod 3) and 2 is not a cubic residue mod p.
    if p % 3 != 1:
        _memo_is_inert[p] = False
        return False
        
    e = (p - 1) // 3
    result = pow(2, e, p) != 1
    _memo_is_inert[p] = result
    return result

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
        if d == 1:
            continue
        if d == n: # Cycle or failure, try a new c
            c += 1
            x, y = 2, 2
            continue
        return d

# --- Main Logic ---

_memo_is_representable = {}
def is_representable(n):
    """
    Recursively checks if n is representable, with memoization.
    """
    if n in _memo_is_representable:
        return _memo_is_representable[n]

    if n == 1:
        return True

    if is_prime(n):
        result = not is_inert(n)
        _memo_is_representable[n] = result
        return result

    # Check for perfect cube, which is always representable (x=c, y=0, z=0)
    # A simple integer root check
    c = round(n**(1/3.0))
    if (c-1)**3 == n or c**3 == n or (c+1)**3 == n:
        _memo_is_representable[n] = True
        return True

    # Find a factor and check its exponent
    factor = pollard_rho(n)
    
    exponent = 0
    temp_n = n
    while temp_n % factor == 0:
        exponent += 1
        temp_n //= factor
        
    if is_inert(factor):
        if exponent % 3 != 0:
            _memo_is_representable[n] = False
            return False
            
    # Recurse on the remaining part of the number
    result = is_representable(temp_n)
    _memo_is_representable[n] = result
    return result

def solve():
    """
    Main function to solve the problem.
    """
    start_n = 10**18
    num_integers_to_check = 10000
    
    count = 0
    
    for k in range(num_integers_to_check + 1):
        n = start_n + k
        if is_representable(n):
            count += 1
            # The prompt asks to "output each number in the final equation"
            # which is interpreted as printing the numbers `n` that satisfy the condition.
            print(n)
    
    print(f"Total count: {count}")

if __name__ == '__main__':
    solve()