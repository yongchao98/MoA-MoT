import numpy as np
from scipy.integrate import quad

# The integral to be solved is: Integral from 0 to infinity of 4 / (1 + 3 * exp(2 * tau^2)) d(tau)
# The numbers in the final equation are:
c1 = 4
c2 = 1
c3 = 3
c4 = 2

print(f"The integral is of the form: integral(c1 / (c2 + c3 * exp(c4 * tau^2)))")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"c4 = {c4}")

# Define the integrand function
def integrand(tau):
    return c1 / (c2 + c3 * np.exp(c4 * tau**2))

# Perform the numerical integration from 0 to infinity
result, error = quad(integrand, 0, np.inf)

print(f"\nThe numerical result of the integral is: {result}")
print(f"The analytical result is pi / sqrt(3).")
print(f"pi / sqrt(3) = {np.pi / np.sqrt(3)}")
