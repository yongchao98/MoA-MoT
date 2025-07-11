import numpy as np
from scipy.integrate import quad
from fractions import Fraction

# Define the simplified integrand
def integrand(x):
    return (np.sin(x) * (2 * np.cos(x) - 1))**50

# Perform numerical integration with high precision
integral_value, error = quad(integrand, 0, np.pi, limit=200)

# The result of such integrals often has a factor of pi.
# We will divide by pi and find the fraction.
value_over_pi = integral_value / np.pi

# Find the closest fraction to this value
frac = Fraction(value_over_pi).limit_denominator(10**10)

numerator = frac.numerator
denominator = frac.denominator

# The requested answer is strictly a fraction. We print the numerator and denominator.
print(f"The result of the integration is approximately {integral_value}.")
print(f"The value divided by pi is approximately {value_over_pi}.")
print(f"This value corresponds to the fraction {numerator}/{denominator}.")
print("Final Equation: ")
print(f"I/pi = {numerator} / {denominator}")

final_answer = str(numerator) + "/" + str(denominator)
print(f"<<<{final_answer}>>>")