from fractions import Fraction

# From the step-by-step derivation, we found the function f(x)
# that satisfies all the given conditions.
# The coefficients are a = -3/8, b = -5/8, and c = 0.
# So, f(x) = x^3 - (3/8)x^2 - (5/8)x.

# We need to compute the exact value of f(3).
x = 3

# The equation for f(3) is:
# f(3) = 3^3 - (3/8) * 3^2 - (5/8) * 3
# f(3) = 27 - 27/8 - 15/8
# f(3) = 27 - 42/8
# f(3) = 27 - 21/4
# f(3) = 108/4 - 21/4
# f(3) = 87/4

# Let's define the numbers in the equation to be printed.
term1 = x**3
term2_num = -27
term2_den = 8
term3_num = -15
term3_den = 8
result_num = 87
result_den = 4

# The problem asks to output each number in the final equation.
# The final numerical equation is 27 - 27/8 - 15/8 = 87/4.
# We print the components of this calculation.
print(f"f(3) = {term1} + ({term2_num}/{term2_den}) + ({term3_num}/{term3_den}) = {result_num}/{result_den}")