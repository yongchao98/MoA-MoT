import numpy as np
from scipy.integrate import quad

# Define the integrand function based on the derived expression.
# The variable 't' in the code corresponds to the variable τ in the integral.
def integrand(t):
    """
    This is the function 1 / (x(0;τ) + y(0;τ) + z(0;τ)) which we need to integrate.
    """
    numerator = 4
    # The values 5, 2, and -1 are from the derived denominator 5*exp(2*τ^2) - 1.
    denominator = 5 * np.exp(2 * t**2) - 1
    return numerator / denominator

# Calculate the definite integral from 0 to infinity.
integral_value, error = quad(integrand, 0, np.inf)

# The problem asks to output the numbers in the final equation.
# The final "equation" is the calculation of the integral I.
# I = integral from 0 to inf of 4 / (5 * exp(2*τ^2) - 1) dτ = result
# The numbers defining this are 4, 5, 2, -1.
print("The integral to be calculated is: I = integral from 0 to inf of_f(τ) dτ")
print("where f(τ) = a / (b * exp(c * τ^2) + d)")
print(f"a = {4}")
print(f"b = {5}")
print(f"c = {2}")
print(f"d = {-1}")

# Print the final numerical result of the integral.
print("\nThe calculated value of the integral is:")
print(integral_value)
