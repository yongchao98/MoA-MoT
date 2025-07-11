from fractions import Fraction

# We have the following system of equations for the expected number of rolls:
# E = 2 + (9/49) * C9
# C9 = 1 + (1/21) * C3
# C3 = 1 + (1/21) * E

# We solve this system step by step.

# First, express C3 in terms of E.
# C3 = 1 + (1/21) * E
c3_const = Fraction(1, 1)
c3_coeff_E = Fraction(1, 21)
print(f"C3 = {c3_const} + {c3_coeff_E} * E")

# Next, substitute the expression for C3 into the equation for C9.
# C9 = 1 + (1/21) * (1 + (1/21) * E)
c9_const = 1 + Fraction(1, 21) * c3_const
c9_coeff_E = Fraction(1, 21) * c3_coeff_E
print(f"C9 = {c9_const} + {c9_coeff_E} * E")

# Finally, substitute the expression for C9 into the equation for E.
# E = 2 + (9/49) * (c9_const + c9_coeff_E * E)
# E = 2 + (9/49)*c9_const + (9/49)*c9_coeff_E * E
# E * (1 - (9/49)*c9_coeff_E) = 2 + (9/49)*c9_const

E_coeff = 1 - Fraction(9, 49) * c9_coeff_E
constant_term = 2 + Fraction(9, 49) * c9_const

# We have the final equation of the form: E * E_coeff = constant_term
# Let's print the numbers in this equation.
print("\nThe final equation is:")
print(f"E * ({E_coeff.numerator}/{E_coeff.denominator}) = {constant_term.numerator}/{constant_term.denominator}")
print("\nNumbers in the final equation before solving for E:")
print(f"Coefficient of E (numerator): {E_coeff.numerator}")
print(f"Coefficient of E (denominator): {E_coeff.denominator}")
print(f"Constant term (numerator): {constant_term.numerator}")
print(f"Constant term (denominator): {constant_term.denominator}")

# Solve for E
E = constant_term / E_coeff

print("\n---")
print(f"The minimal expected value of rolls is a fraction: {E}")
print(f"Numerator: {E.numerator}")
print(f"Denominator: {E.denominator}")