import numpy as np
from scipy.integrate import quad

def p(x):
    """
    Defines the piecewise function p(x).
    """
    if 0 <= x <= 3:
        return (2 * x**3) / 8
    elif 3 < x <= 5:
        return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))
    else:
        return 0

# Since the integration range is from 0 to 4, we need to split the integral
# at x=3 where the function definition changes.

# Part 1: Integral from 0 to 3
integral_part1, error1 = quad(lambda x: (2 * x**3) / 8, 0, 3)

# Part 2: Integral from 3 to 4
integral_part2, error2 = quad(lambda x: (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x)), 3, 4)

# Total integral is the sum of the two parts
total_integral = integral_part1 + integral_part2

print(f"The integral is split into two parts at x=3:")
print(f"First part (from 0 to 3): {integral_part1}")
print(f"Second part (from 3 to 4): {integral_part2}")
print("The final equation is:")
print(f"{integral_part1} + ({integral_part2}) = {total_integral}")
print(f"The total integral from x=0 to x=4 is {total_integral}")
