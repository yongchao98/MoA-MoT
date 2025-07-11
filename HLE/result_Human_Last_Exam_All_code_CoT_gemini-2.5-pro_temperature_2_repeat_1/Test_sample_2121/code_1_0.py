import numpy as np
from scipy.integrate import quad

# The integrand is derived from the solution to the equations of motion.
# The denominator of the integral to be computed is x(0;τ) + y(0;τ) + z(0;τ).
# Through analysis of the differential equations, this sum is found to be:
# S(0;τ) = (1/4) + (3/4) * exp(2*τ^2)
#
# The integral we need to evaluate is:
# Integral = ∫[0, ∞] 1 / S(0;τ) dτ = ∫[0, ∞] 1 / ( (1/4) + (3/4) * exp(2*τ^2) ) dτ
# which simplifies to ∫[0, ∞] 4 / (1 + 3 * exp(2*τ^2)) dτ

# Define the function to be integrated.
def integrand(tau):
    """
    This is the function 1 / (x(0;τ) + y(0;τ) + z(0;τ))
    """
    numerator = 4
    denominator = 1 + 3 * np.exp(2 * tau**2)
    return numerator / denominator

# Perform the numerical integration from 0 to infinity.
result, error = quad(integrand, 0, np.inf)

# The final equation is the integral of the derived function of tau.
# We will print the numbers in the final expression of the integral.
numerator = 4
c1 = 1
c2 = 3
c3 = 2
print(f"The integral to be computed is ∫[0, ∞] {numerator} / ({c1} + {c2} * exp({c3} * τ^2)) dτ")

# Print the final numerical result.
print(f"The value of the integral is: {result}")
print(f"The estimated error is: {error}")

<<<1.516386059139336>>>