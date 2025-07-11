import math

def prime_factorization(n):
    """
    Finds the prime factorization of a number n.
    Returns a dictionary of {prime: exponent}.
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

def phi_from_factors(n, factors):
    """
    Calculates Euler's totient function phi(n) using its prime factorization.
    """
    if n == 1:
        return 1
    result = n
    for p in factors:
        result = result // p * (p - 1)
    return result

# Step 1: Define n
n = math.factorial(10)

# Step 2: Find the prime factorization of n
# The only prime factors of 10! are primes <= 10, which are 2, 3, 5, 7.
factors_of_n = prime_factorization(n)

# Step 3: Calculate phi(n)
phi_n = phi_from_factors(n, factors_of_n)

# Step 4: Calculate the number of manifolds
# For n > 2, the number of manifolds is phi(n) / 2
num_manifolds = phi_n // 2

# Final Output
print(f"The number of closed orientable 3-manifolds with a fundamental group of a given finite order n is found by counting specific quotient spaces of the 3-sphere.")
print(f"Let n be the order of the fundamental group.")
print(f"n = 10! = {n}")
print(f"Under the assumption that the only such fundamental group is the cyclic group C_n, the number of manifolds is given by phi(n) / 2.")
print(f"The final calculation is: phi({n}) / 2")
print(f"phi({n}) = {phi_n}")
print(f"Number of manifolds = {phi_n} / 2 = {num_manifolds}")