import numpy as np
from scipy.integrate import quad
from fractions import Fraction

def integrand(x):
    """
    This is the simplified integrand function derived from the original problem.
    Original: max(|2*sin(x)|, |2*cos(2*x) - 1|)**50 * min(|sin(2*x)|, |cos(3*x)|)**50
    Simplified: (|sin(2*x) * (2*cos(2*x) - 1)|)**50
    """
    term1 = np.sin(2 * x)
    term2 = 2 * np.cos(2 * x) - 1
    
    # The even power of 50 makes the absolute value redundant for the final value,
    # but it's good practice to keep it for clarity with the derivation.
    return (np.abs(term1 * term2))**50

# Define the integration limits
lower_limit = 0
upper_limit = np.pi

# Perform the numerical integration
result, error = quad(integrand, lower_limit, upper_limit)

# The result of the integration is a floating-point number.
# We will convert it to the closest fraction as requested.
# The limit_denominator() method is used to find a rational approximation.
fraction_result = Fraction(result).limit_denominator()

# As requested, printing the numbers in the final equation (i.e., the fraction).
print(f"The result of the integration is: {fraction_result}")
print(f"Numerator: {fraction_result.numerator}")
print(f"Denominator: {fraction_result.denominator}")
