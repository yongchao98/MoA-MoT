from fractions import Fraction

# We derived the following system of linear equations for the expected values:
# 1) E_1 = 2 + (9/49) * E_9
# 2) E_9 = 1 + (1/21) * E_3
# 3) E_3 = 1 + (1/21) * E_1
# We want to solve for E_1.

# Let's use substitution. First, express E_3 in terms of E_1 from equation (3).
# E_3 = 1 + (1/21) * E_1
e3_coeff_e1 = Fraction(1, 21)
e3_const = Fraction(1, 1)

# Now, substitute this expression for E_3 into equation (2) to get E_9 in terms of E_1.
# E_9 = 1 + (1/21) * (1 + (1/21) * E_1)
# E_9 = 1 + 1/21 + (1/21)^2 * E_1
# E_9 = 22/21 + (1/441) * E_1
e9_coeff_e1 = Fraction(1, 21) * e3_coeff_e1
e9_const = Fraction(1, 1) + Fraction(1, 21) * e3_const

# Finally, substitute the expression for E_9 into equation (1).
# E_1 = 2 + (9/49) * E_9
# E_1 = 2 + (9/49) * (e9_const + e9_coeff_e1 * E_1)
# E_1 = 2 + (9/49) * e9_const + (9/49) * e9_coeff_e1 * E_1

# Let's calculate the terms.
const_term = Fraction(2, 1) + Fraction(9, 49) * e9_const
e1_coeff = Fraction(9, 49) * e9_coeff_e1

# E_1 = const_term + e1_coeff * E_1
# E_1 * (1 - e1_coeff) = const_term
# E_1 = const_term / (1 - e1_coeff)

e1 = const_term / (Fraction(1, 1) - e1_coeff)

# Print the final result
# The variable 'e1' holds the final answer as a Fraction object.
# We print its components to show the simplified fraction.
print("The minimal expected number of rolls E is given by the equation:")
print(f"E = {e1.numerator} / {e1.denominator}")