import numpy as np
from scipy.integrate import quad

# Define the integrand function
def integrand(tau):
    return 4 / (1 + 3 * np.exp(2 * tau**2))

# Calculate the definite integral from 0 to infinity
integral_value, error = quad(integrand, 0, np.inf)

# The integral can be written in the form:
# integral from 0 to infinity of A / (B + C * exp(D * t**2)) dt
A = 4
B = 1
C = 3
D = 2

# Output the numbers in the final equation as requested
print(f"The integral is of the form: integral(A / (B + C * exp(D * t**2)))")
print(f"The parameters are: A = {A}, B = {B}, C = {C}, D = {D}")

# Print the final numerical result
print("\nThe value of the time-averaged integral is:")
print(integral_value)
