# This script constructs and prints the derived formula for P(n).
# The formula is composed of two terms, corresponding to the n^-2 and n^-3 corrections.
# L represents ln(n).

# Coefficients for the n^-2 term
num1_coeff_L2 = 3
num1_coeff_L1 = -2
num1_coeff_L0 = 2
den1_coeff_n2 = 24

# Coefficients for the n^-3 term
num2_coeff_L3 = 1
num2_coeff_L2 = -2
num2_coeff_L1 = 2
den2_coeff_n3 = 48

# Construct the formula as a string for display.
# Using '^' for exponentiation to represent the mathematical formula.
term1_numerator = f"({num1_coeff_L2}*L^2 - {abs(num1_coeff_L1)}*L + {num1_coeff_L0})"
term1_denominator = f"({den1_coeff_n2}*n^2)"
term1 = f"{term1_numerator}/{term1_denominator}"

term2_numerator = f"(L^3 - {abs(num2_coeff_L2)}*L^2 + {num2_coeff_L1}*L)"
term2_denominator = f"({den2_coeff_n3}*n^3)"
term2 = f"{term2_numerator}/{term2_denominator}"

# The full formula for P(n) is the sum of these two terms.
final_formula = f"{term1} + {term2}"

print(final_formula)