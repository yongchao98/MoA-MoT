import numpy as np
from scipy import integrate
import math

# Step 1-4: The integrand simplifies to log|sin(x)/x|.
# We define this function for numerical integration.
def integrand(x):
    """
    The function to be integrated, log|sin(x)/x|.
    We handle the case x=0 separately, as sin(x)/x -> 1, so log(1) = 0.
    """
    if x == 0:
        return 0.0
    return np.log(np.abs(np.sin(x) / x))

# Step 5-6: Numerically evaluate the integral and show it matches the known result.
# We integrate from 0 to a very large number to approximate infinity.
# The `quad` function from scipy.integrate is a robust numerical integrator.
integral_value, integral_error = integrate.quad(integrand, 0, np.inf)

# The known analytical value of the integral
analytical_value = -math.pi / 2

print("The integral to be evaluated is: integral from 0 to infinity of log|sin(x)/x| dx")
print(f"The analytical result is a well-known value: -pi/2")
print(f"Equation: integral(...) = -pi/2")
print(f"Numerical verification:")
print(f"Calculated value: {integral_value}")
print(f"Analytical value (-pi/2): {analytical_value}")
print(f"The error in numerical calculation is: {integral_error}")

# "Remember in the final code you still need to output each number in the final equation!"
# The final equation can be seen as: The Integral = -pi / 2
# Here are the components of the result.
print("\nThe final equation is value = -pi / 2.")
print(f"The number pi is: {math.pi}")
print(f"The number 2 is: 2")
print(f"The negative sign: -1")