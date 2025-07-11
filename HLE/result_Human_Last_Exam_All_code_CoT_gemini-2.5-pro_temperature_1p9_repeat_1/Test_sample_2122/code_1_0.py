import numpy as np
from scipy.integrate import quad

# The final temperature T is given by the formula:
# T = -u + exp(u^2) / (e + integral)
# where u = 2 and integral is the definite integral of exp(x^2) from 1 to u.

u = 2

# Calculate the value of the definite integral
integral_value, integral_error = quad(lambda x: np.exp(x**2), 1, u)

# Values of the constants
e = np.e
exp_u_squared = np.exp(u**2)

# Calculate the denominator of the fraction
denominator = e + integral_value

# Calculate the final temperature
temperature = -u + exp_u_squared / denominator

# Print out the numbers used in the final equation step by step
print(f"The calculation is for u = {u}")
print(f"The equation for the temperature is: T(u) = -u + exp(u^2) / (e + integral(exp(x^2) from 1 to u))")
print(f"Value of the integral from 1 to {u}: {integral_value}")
print(f"Value of e: {e}")
print(f"Value of exp({u}^2): {exp_u_squared}")
print(f"The term (e + integral): {denominator}")
print(f"The fractional part exp(u^2) / (e + integral): {exp_u_squared / denominator}")
print(f"Final temperature T({u}): {temperature}")

<<<1.0819803154162812>>>