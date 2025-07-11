import math

def get_prime_factorization(n):
    """
    Returns a list of prime factors of n.
    """
    factors = []
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
       factors.append(temp)
    return factors

def gcd(a, b):
    """
    Computes the greatest common divisor of a and b.
    """
    return math.gcd(a, b)

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.")
    print("-" * 60)

    # Step 1: Prime factorization of d
    primes = get_prime_factorization(d)
    print(f"Step 1: The prime factorization of d = {d} is: {' * '.join(map(str, primes))}")
    print("Since d is square-free, a character is primitive if its components are non-trivial.")
    print("-" * 60)

    # Step 2: Calculate p-1 for each prime factor
    p_minus_1 = [p - 1 for p in primes]
    print("Step 2: For each prime factor p, calculate p-1:")
    for i in range(len(primes)):
        print(f"  For p = {primes[i]}, p-1 = {p_minus_1[i]}")
    print("-" * 60)

    # Step 3: Use inclusion-exclusion principle
    print("Step 3: Use inclusion-exclusion to find characters of order exactly 6.")
    print("The number is N = |A_6| - |A_3| - |A_2| + |A_1|, where |A_k| is the")
    print("number of primitive characters with order dividing k.")
    print("-" * 60)

    # Step 4: Calculate |A_k| for divisors of g=6
    divisors_of_g = [1, 2, 3, 6]
    A = {}
    
    print("Step 4: Calculate |A_k| for each divisor k of 6:")
    for k in sorted(divisors_of_g):
        terms = []
        prod = 1
        print(f"\nCalculating |A_{k}| = PRODUCT(gcd(k, p-1) - 1):")
        for i in range(len(primes)):
            p = primes[i]
            pm1 = p_minus_1[i]
            val = gcd(k, pm1) - 1
            terms.append(val)
            prod *= val
            print(f"  p={p:<2}: gcd({k}, {pm1}) - 1 = {val}")
        A[k] = prod
        print(f"Result: |A_{k}| = {' * '.join(map(str, terms))} = {prod}")

    print("-" * 60)

    # Step 5: Calculate the final result
    print("Step 5: Apply the inclusion-exclusion formula:")
    result = A[6] - A[3] - A[2] + A[1]
    
    print(f"Number of characters of order 6 = |A_6| - |A_3| - |A_2| + |A_1|")
    print(f"                               = {A[6]} - {A[3]} - {A[2]} + {A[1]}")
    print(f"                               = {result}")
    
    print("-" * 60)
    print(f"The final number of primitive Dirichlet characters is {result}.")

solve()