import math

def get_prime_factorization(n):
    """
    Computes the prime factorization of a given integer n.
    Returns a dictionary where keys are prime factors and values are their exponents.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
       factors[temp] = factors.get(temp, 0) + 1
    return factors

def tau_from_factors(factors):
    """
    Computes the number of divisors (tau function) from prime factorization.
    """
    num_divisors = 1
    for p in factors:
        num_divisors *= (factors[p] + 1)
    return num_divisors

def phi_from_factors(n, factors):
    """
    Computes Euler's totient function (phi function) from prime factorization.
    """
    result = n
    for p in factors:
        result = result * (p - 1) // p
    return result

def solve():
    """
    Calculates the number of closed orientable 3-manifolds with a cyclic
    fundamental group of cardinality 10!, based on the formula for lens spaces.
    """
    n = math.factorial(10)
    
    # Get prime factorization of 10!
    # Primes <= 10 are 2, 3, 5, 7
    factors = {2: 0, 3: 0, 5: 0, 7: 0}
    for i in range(2, 11):
        num_factors = get_prime_factorization(i)
        for p, exp in num_factors.items():
            factors[p] = factors.get(p, 0) + exp

    # Calculate tau(10!) and phi(10!)
    tau_n = tau_from_factors(factors)
    phi_n = phi_from_factors(n, factors)
    
    # For n > 2, the number of non-homeomorphic lens spaces L(n,q)
    # is given by (tau(n) + phi(n)) / 2
    if n > 2:
        num_manifolds = (tau_n + phi_n) // 2
    else: # n=1 or n=2 are special cases
        num_manifolds = 1

    print("This solution assumes we are counting manifolds for the simplest case where the fundamental group is the cyclic group C(10!).")
    print(f"The order of the fundamental group is n = 10! = {n}.")
    print(f"The number of divisors is tau(n) = {tau_n}.")
    print(f"Euler's totient function is phi(n) = {phi_n}.")
    print("\nThe number of such manifolds (lens spaces) is calculated by the formula (phi(n) + tau(n)) / 2:")
    print(f"({phi_n} + {tau_n}) / 2 = {num_manifolds}")

solve()
