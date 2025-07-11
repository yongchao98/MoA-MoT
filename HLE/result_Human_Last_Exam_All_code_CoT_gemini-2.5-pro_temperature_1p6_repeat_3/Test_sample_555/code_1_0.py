from fractions import Fraction

# We want to find E_1, the expected number of rolls starting with no information.
# The system of equations for E_1, E_3, E_7, E_9 is:
# E1 = 1 + E7
# E7 = 1 + (9/49) * E9
# E9 = 1 + (3/63) * E3 = 1 + (1/21) * E3
# E3 = 1 + (1/21) * E1

# We can solve this system by substitution.
# All calculations will be done using the Fraction class for exact results.

# E3 = 1 + (1/21) * E1
E3_const = Fraction(1)
E3_coeff_E1 = Fraction(1, 21)
print(f"Equation for E3: E3 = {E3_const} + {E3_coeff_E1} * E1")

# E9 = 1 + (1/21) * E3
# Substitute E3 into the equation for E9
E9_const = 1 + Fraction(1, 21) * E3_const
E9_coeff_E1 = Fraction(1, 21) * E3_coeff_E1
print(f"Equation for E9: E9 = {E9_const} + {E9_coeff_E1} * E1")

# E7 = 1 + (9/49) * E9
# Substitute E9 into the equation for E7
E7_const = 1 + Fraction(9, 49) * E9_const
E7_coeff_E1 = Fraction(9, 49) * E9_coeff_E1
print(f"Equation for E7: E7 = {E7_const} + {E7_coeff_E1} * E1")

# E1 = 1 + E7
# Substitute E7 into the equation for E1
E1_const_from_E7 = 1 + E7_const
E1_coeff_from_E7 = E7_coeff_E1
print(f"Equation for E1: E1 = {E1_const_from_E7} + {E1_coeff_from_E7} * E1")

# Now, solve for E1:
# E1 = E1_const_from_E7 + E1_coeff_from_E7 * E1
# E1 * (1 - E1_coeff_from_E7) = E1_const_from_E7
# E1 = E1_const_from_E7 / (1 - E1_coeff_from_E7)

lhs_coeff = 1 - E1_coeff_from_E7
rhs_const = E1_const_from_E7

final_E1 = rhs_const / lhs_coeff

print("\nThe final equation to solve for E1 is:")
print(f"E1 * (1 - {E1_coeff_from_E7.numerator}/{E1_coeff_from_E7.denominator}) = {E1_const_from_E7.numerator}/{E1_const_from_E7.denominator}")
print(f"E1 * ({lhs_coeff.numerator}/{lhs_coeff.denominator}) = {rhs_const.numerator}/{rhs_const.denominator}")

print("\nSolving for E1:")
print(f"E1 = ({rhs_const.numerator}/{rhs_const.denominator}) / ({lhs_coeff.numerator}/{lhs_coeff.denominator})")
print(f"E1 = {final_E1.numerator}/{final_E1.denominator}")
print(f"\nThe minimal expected number of rolls is {final_E1.numerator}/{final_E1.denominator}.")
