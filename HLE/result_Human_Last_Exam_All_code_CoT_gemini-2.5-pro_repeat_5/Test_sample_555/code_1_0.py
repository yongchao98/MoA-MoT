from fractions import Fraction

# We have the following system of equations for the expected number of rolls E(n):
# E(1) = 1 + E(7)
# E(7) = 1 + (9/49) * E(9)
# E(9) = 1 + (1/21) * E(3)
# E(3) = 1 + (1/21) * E(1)
#
# We want to solve for E(1). Let's use substitution.
# Let E1, E3, E7, E9 be the expected values.
#
# From E(3) = 1 + (1/21) * E(1), substitute into the equation for E(9):
# E(9) = 1 + Fraction(1, 21) * (1 + Fraction(1, 21) * E(1))
# E(9) = 1 + Fraction(1, 21) + Fraction(1, 441) * E(1)
# E(9) = Fraction(22, 21) + Fraction(1, 441) * E(1)
#
# Substitute this into the equation for E(7):
# E(7) = 1 + Fraction(9, 49) * (Fraction(22, 21) + Fraction(1, 441) * E(1))
# E(7) = 1 + Fraction(9 * 22, 49 * 21) + Fraction(9 * 1, 49 * 441) * E(1)
# E(7) = 1 + Fraction(198, 1029) + Fraction(9, 21609) * E(1)
# Simplify fractions: 198/1029 = 66/343; 9/21609 = 1/2401
# E(7) = 1 + Fraction(66, 343) + Fraction(1, 2401) * E(1)
#
# Finally, substitute this into the equation for E(1):
# E(1) = 1 + E(7)
# E(1) = 1 + (1 + Fraction(66, 343) + Fraction(1, 2401) * E(1))
# E(1) = 2 + Fraction(66, 343) + Fraction(1, 2401) * E(1)
#
# Now, we solve this final equation for E(1).
# E(1) - Fraction(1, 2401) * E(1) = 2 + Fraction(66, 343)
# E(1) * (1 - Fraction(1, 2401)) = 2 + Fraction(66, 343)

# Represent the equation with Fraction objects
lhs_coeff = Fraction(1, 1) - Fraction(1, 2401)
rhs_value = Fraction(2, 1) + Fraction(66, 343)

# Solve for E1
# E1 * lhs_coeff = rhs_value  =>  E1 = rhs_value / lhs_coeff
E1 = rhs_value / lhs_coeff

# The final fraction is simplified by the Fraction class automatically.
numerator = E1.numerator
denominator = E1.denominator

# Output the result as an equation
print("The minimal expected value of rolls, E, is given by the solution to a system of linear equations.")
print("The final simplified equation for E is:")
print(f"E = {numerator}/{denominator}")
