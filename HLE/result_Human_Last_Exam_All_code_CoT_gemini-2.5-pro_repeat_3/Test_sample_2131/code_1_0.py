import math

# The derived analytical expression for the deflection y(0) is (3/4) * (9/4)^(3/5).
# Here, we define the numbers in this final equation and compute the result.

# Numerator and denominator of the coefficient
coeff_num = 3
coeff_den = 4

# Numerator and denominator of the base
base_num = 9
base_den = 4

# Numerator and denominator of the exponent
exp_num = 3
exp_den = 5

# Calculate the terms
coefficient = coeff_num / coeff_den
base = base_num / base_den
exponent = exp_num / exp_den

# Calculate the final value of y(0)
y_at_0 = coefficient * (base ** exponent)

# Print the components of the equation and the final result
print(f"The equation for y(0) is ({coeff_num}/{coeff_den}) * ({base_num}/{base_den})^({exp_num}/{exp_den})")
print(f"The membrane's deflection at x = 0, y(0), is: {y_at_0}")