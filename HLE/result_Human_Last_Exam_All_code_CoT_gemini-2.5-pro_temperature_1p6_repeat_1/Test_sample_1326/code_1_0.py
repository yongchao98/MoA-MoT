import numpy as np
from scipy.integrate import quad

# Define the first part of the function for 0 <= x <= 3
def p1(x):
    return (x**3) / 4

# Define the second part of the function for 3 <= x <= 5
def p2(x):
    # Ensure cosine is not -1 to avoid division by zero, though not an issue in the integration range [3, 4]
    if np.cos(x) == -1.0:
        return np.inf 
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral of the first part from 0 to 3
# The result of quad is a tuple (integral_value, error_estimate)
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral of the second part from 3 to 4
integral2, error2 = quad(p2, 3, 4)

# Calculate the total integral
total_integral = integral1 + integral2

# Print the final equation with the calculated values
print(f"The integral from 0 to 4 is the sum of the integral from 0 to 3 and the integral from 3 to 4.")
print(f"The equation is: {integral1} + ({integral2}) = {total_integral}")
