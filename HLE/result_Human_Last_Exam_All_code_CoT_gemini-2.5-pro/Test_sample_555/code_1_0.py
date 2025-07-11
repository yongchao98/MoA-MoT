from fractions import Fraction

# We have the following system of equations:
# 1) E = 2 + (9/49) * E_9
# 2) E_9 = 1 + (1/21) * E_3
# 3) E_3 = 1 + (1/21) * E
# We want to solve for E.

# Let's express the coefficients as Fraction objects for exact arithmetic.
c1 = Fraction(9, 49)
c2 = Fraction(1, 21)
c3 = Fraction(1, 21)

print("The system of equations for the expected values is:")
print(f"E = 2 + ({c1.numerator}/{c1.denominator}) * E_9")
print(f"E_9 = 1 + ({c2.numerator}/{c2.denominator}) * E_3")
print(f"E_3 = 1 + ({c3.numerator}/{c3.denominator}) * E")
print("-" * 20)

# Step 1: Substitute E_3 from equation (3) into equation (2).
# E_9 = 1 + (1/21) * (1 + (1/21) * E)
# E_9 = 1 + 1/21 + (1/441) * E
# E_9 = 22/21 + (1/441) * E
E9_in_terms_of_E_const = Fraction(1) + c2 * Fraction(1)
E9_in_terms_of_E_coeff = c2 * c3
print("Solving for E_9 in terms of E:")
print(f"E_9 = {E9_in_terms_of_E_const.numerator}/{E9_in_terms_of_E_const.denominator} + ({E9_in_terms_of_E_coeff.numerator}/{E9_in_terms_of_E_coeff.denominator}) * E")
print("-" * 20)

# Step 2: Substitute the expression for E_9 into equation (1).
# E = 2 + (9/49) * (22/21 + (1/441) * E)
# E = 2 + (9/49)*(22/21) + (9/49)*(1/441) * E
# E = 2 + 198/1029 + 9/21609 * E
# E * (1 - 9/21609) = 2 + 198/1029
# E * (2400/2401) = 752/343
# E = (752/343) * (2401/2400)

# Calculate the constant term on the right side
rhs_const = Fraction(2) + c1 * E9_in_terms_of_E_const
# Calculate the coefficient of E on the right side
rhs_E_coeff = c1 * E9_in_terms_of_E_coeff

# Rearrange to have E terms on one side: E * (1 - rhs_E_coeff) = rhs_const
E_coeff_final = Fraction(1) - rhs_E_coeff

print("Solving for E:")
# Print the penultimate equation
print(f"E * ({E_coeff_final.numerator}/{E_coeff_final.denominator}) = {rhs_const.numerator}/{rhs_const.denominator}")

# Solve for E
E = rhs_const / E_coeff_final

print("-" * 20)
print("The final calculation is:")
# Multiplying by the reciprocal
reciprocal = Fraction(E_coeff_final.denominator, E_coeff_final.numerator)
print(f"E = ({rhs_const.numerator}/{rhs_const.denominator}) * ({reciprocal.numerator}/{reciprocal.denominator})")
print(f"E = {E.numerator}/{E.denominator}")
print("-" * 20)
print("The minimal expected number of rolls is:")
print(f"{E.numerator}/{E.denominator}")
