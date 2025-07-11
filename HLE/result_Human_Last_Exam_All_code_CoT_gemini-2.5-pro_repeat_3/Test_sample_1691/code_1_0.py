import math

# This script develops and presents an analytical approximation for the integral
# I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8) dx
# for the small epsilon regime.

# 1. Identify parameters from the dominant term near x=0.
# The denominator's polynomial is g(x) = 9.0*x**5.0 + 5.0*x**6.0 + 9.0*x**8.0.
# For small x, g(x) is approximated by its leading term, c0*x**p.
c0 = 9.0
p = 5.0

# 2. The approximation for I(epsilon) has the form C * epsilon^alpha.
# Calculate the exponent alpha.
alpha = -(p - 1.0) / p

# 3. Calculate the coefficient C using the derived formula.
# C = (c0**(-1/p) / p) * (pi / sin(pi/p))
C = (c0**(-1.0/p) / p) * (math.pi / math.sin(math.pi/p))

# 4. Print the results in a structured way.
print("The analytical formula for the integral I(epsilon) for small epsilon is developed as follows:")
print("-" * 70)
print("1. The integral is approximated by analyzing the integrand's behavior near x=0.")
print("   The function g(x) in the denominator is approximated by its dominant term:")
print(f"   g(x) ≈ c0 * x**p, with c0 = {c0}, p = {p}")
print("")

print("2. This leads to an approximate formula of the form: I(epsilon) ≈ C * epsilon**alpha")
print("")

print("3. The exponent 'alpha' is determined by p:")
print(f"   alpha = - (p - 1) / p = - ({p} - 1.0) / {p} = {alpha}")
print("")

print("4. The coefficient 'C' is calculated from c0 and p:")
print("   The formula for C is: (c0**(-1/p) / p) * (pi / sin(pi/p))")
# Outputting each number in the final equation for C
c0_term = c0**(-1.0/p)
p_term = p
pi_term = math.pi
sin_term = math.sin(math.pi/p)
print(f"   - c0**(-1/p) = {c0}**(-1.0/{p}) = {c0_term}")
print(f"   - p          = {p_term}")
print(f"   - pi         = {pi_term}")
print(f"   - sin(pi/p)  = sin(pi/{p}) = {sin_term}")
print(f"   Combining these values: C = ({c0_term} / {p_term}) * ({pi_term} / {sin_term})")
print(f"   So, C = {C}")
print("-" * 70)

print("The final analytical approximation is:")
print(f"I(epsilon) ≈ {C:.5f} * epsilon**({alpha})")
print("-" * 70)
