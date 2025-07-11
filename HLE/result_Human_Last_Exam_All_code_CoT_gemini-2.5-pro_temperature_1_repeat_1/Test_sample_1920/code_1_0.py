import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
       factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def get_divisors(n):
    """Returns a sorted list of all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def euler_phi(n):
    """Computes Euler's totient function phi(n)."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

def mobius_mu(n):
    """Computes the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Find the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.")
    
    d_factors = get_prime_factorization(d)
    primes = list(d_factors.keys())
    print(f"The conductor d = {d} is square-free, with prime factors: {primes}.")
    print("A character chi modulo d is primitive iff its p-components are non-trivial for all p|d.")
    
    print("\nThe number of primitive characters of order exactly g, N(g), is found using Mobius inversion:")
    print("N(g) = sum_{k|g} mu(g/k) * C(k)")
    print("where C(k) is the number of primitive characters with order dividing k.")
    print(f"For g=6, the formula is: N(6) = C(6) - C(3) - C(2) + C(1).\n")
    
    c_vals = {}
    g_divisors = get_divisors(g)

    for k in sorted(g_divisors, reverse=True):
        print(f"Calculating C({k}):")
        
        num_chars_per_prime = []
        for p in primes:
            p_minus_1 = p - 1
            limit = gcd(k, p_minus_1)
            limit_divisors = get_divisors(limit)
            
            num_chars_p = sum(euler_phi(j) for j in limit_divisors if j > 1)
            num_chars_per_prime.append(num_chars_p)

        c_k = math.prod(num_chars_per_prime)
        c_vals[k] = c_k

        print(f"  For each prime factor p of d, the number of non-trivial characters mod p with order dividing {k} is sum(phi(j)) for j|gcd({k}, p-1) and j>1.")
        # Show calculation for the first prime as an example
        p, p_minus_1 = primes[0], primes[0]-1
        limit = gcd(k, p_minus_1)
        limit_divisors_gt1 = [j for j in get_divisors(limit) if j > 1]
        if limit_divisors_gt1:
             phi_calcs = ' + '.join([f"phi({j})" for j in limit_divisors_gt1])
             phi_vals = ' + '.join([str(euler_phi(j)) for j in limit_divisors_gt1])
             print(f"  For p={p}, this is {phi_calcs} = {phi_vals} = {num_chars_per_prime[0]}.")
        else:
            print(f"  For p={p}, there are no non-trivial characters with order dividing {k}, so the count is 0.")

        print(f"  The counts for all primes {primes} are {num_chars_per_prime}.")
        print(f"  Therefore, C({k}) = {' * '.join(map(str, num_chars_per_prime))} = {c_k}")

    c1, c2, c3, c6 = c_vals.get(1, 0), c_vals.get(2, 0), c_vals.get(3, 0), c_vals.get(6, 0)
    total_count = c6 - c3 - c2 + c1
    
    print("\nFinally, we plug these values into the formula:")
    print(f"Number of characters = C(6) - C(3) - C(2) + C(1)")
    print(f"                     = {c6} - {c3} - {c2} + {c1}")
    print(f"                     = {total_count}")
    
    print(f"<<<{total_count}>>>")

solve()