import numpy as np
from scipy.integrate import quad

# The final temperature T is given by the formula:
# T = -2 + e^4 / (e + I)
# where I is the definite integral of e^(s^2) from 1 to 2.

# We will now calculate the components of this equation.

# Define the function to integrate
def integrand(s):
    return np.exp(s**2)

# Calculate the values of e and e^4
e_val = np.exp(1)
e4_val = np.exp(4)

# Calculate the definite integral I
integral_val, integral_err = quad(integrand, 1, 2)

# Now, substitute these values back into the equation for T
temperature = -2 + e4_val / (e_val + integral_val)

# Print the components of the equation
print("The final temperature is calculated using the equation: T = A + B / (C + D)")
print(f"Value of A (which is -2): {-2}")
print(f"Value of B (which is e^4): {e4_val}")
print(f"Value of C (which is e): {e_val}")
print(f"Value of D (the integral): {integral_val}")

# Print the final result
print("\nFinal calculated temperature:")
print(temperature)
