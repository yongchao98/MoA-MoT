import math
from fractions import Fraction

# The dimension of the Frostman measure is given.
alpha_num = 8
alpha_den = 5
alpha = Fraction(alpha_num, alpha_den)

# The problem asks for the smallest possible c such that the L^2 norm of the
# Fourier transform of the measure on a circle of radius r is O(r^(c+epsilon)).
# This is a classic result in Fourier restriction theory. The sharp exponent c
# is given by the formula c = -(2 - alpha) / 2, where alpha is the dimension
# of the measure in R^2.

# We will now calculate this value for alpha = 8/5.
# The calculation is c = -(2 - 8/5) / 2

# Step 1: Calculate the term in the parenthesis (2 - alpha)
expr_in_paren = 2 - alpha
# expr_in_paren = 2/1 - 8/5 = 10/5 - 8/5 = 2/5

# Step 2: Negate it and divide by 2
c_frac = -expr_in_paren / 2
# c_frac = -(2/5) / 2 = -1/5

# Now, we print the calculation step-by-step, showing each number.
print(f"The dimension of the measure is alpha = {alpha.numerator}/{alpha.denominator}.")
print("The smallest possible value for c is given by the formula from Fourier restriction theory:")
print("c = -(2 - alpha) / 2")
print("\nSubstituting the value of alpha:")
print(f"c = -(2 - {alpha.numerator}/{alpha.denominator}) / 2")
print(f"c = -({2 * alpha.denominator}/{alpha.denominator} - {alpha.numerator}/{alpha.denominator}) / 2")
print(f"c = -({expr_in_paren.numerator}/{expr_in_paren.denominator}) / 2")
print(f"c = {c_frac.numerator}/{c_frac.denominator}")

# Final answer as a decimal
c_decimal = float(c_frac)
print(f"\nThe value of c is {c_decimal}.")
