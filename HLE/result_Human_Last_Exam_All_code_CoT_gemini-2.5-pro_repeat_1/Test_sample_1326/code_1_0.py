import numpy as np
from scipy.integrate import quad

# Define the two parts of the function p(x)
def part1(x):
    """Function p(x) for 0 <= x <= 3"""
    return (2 * x**3) / 8

def part2(x):
    """Function p(x) for 3 <= x <= 5"""
    # Ensure we are using floating point division
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral for the first interval [0, 3]
# The quad function returns the integral result and an estimated error. We only need the result.
integral1, error1 = quad(part1, 0, 3)

# Calculate the integral for the second interval [3, 4]
integral2, error2 = quad(part2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the final equation with the calculated values
print("The integral of p(x) from 0 to 4 is the sum of the integral from 0 to 3 and the integral from 3 to 4.")
print(f"Integral from 0 to 3: {integral1}")
print(f"Integral from 3 to 4: {integral2}")
print("\nFinal calculation:")
print(f"{integral1} + {integral2} = {total_integral}")
