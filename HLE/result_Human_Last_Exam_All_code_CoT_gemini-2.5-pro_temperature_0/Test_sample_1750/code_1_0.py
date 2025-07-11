import numpy as np
from scipy.integrate import quad

# The simplified integrand is (2 * sin(x) * cos(3*x))^50
# which is also (sin(4*x) - sin(2*x))^50
def integrand(x):
    return (2 * np.sin(x) * np.cos(3 * x))**50

# Integrate from 0 to pi
integral_value, error = quad(integrand, 0, np.pi)

# The final equation is I = integral from 0 to pi of (2*sin(x)*cos(3*x))^50 dx
# The numbers in the equation are 0, pi, 2, 3, 50
print("The simplified integral to be solved is:")
print("I = integral from 0 to pi of (2*sin(x)*cos(3*x))^50 dx")
print("\nThis integral evaluates to a number involving pi, so it cannot be expressed as a simple fraction.")
print("The numerical value is approximately:")
print(integral_value)
