import numpy as np
from scipy.integrate import quad

# Define the first part of the function p(x) for 0 <= x <= 3
def p1(x):
    return (2 * x**3) / 8

# Define the second part of the function p(x) for 3 < x <= 5
def p2(x):
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral of the first part from 0 to 3
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral of the second part from 3 to 4
integral2, error2 = quad(p2, 3, 4)

# Calculate the total integral
total_integral = integral1 + integral2

# Print the final equation with the calculated values
print(f"Integral from 0 to 4 = (Integral from 0 to 3) + (Integral from 3 to 4)")
print(f"Integral = {integral1:.4f} + ({integral2:.4f})")
print(f"Total Integral = {total_integral:.4f}")
