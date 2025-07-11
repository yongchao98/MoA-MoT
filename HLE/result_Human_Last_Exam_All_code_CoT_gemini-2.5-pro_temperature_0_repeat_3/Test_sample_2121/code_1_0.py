import numpy as np
from scipy.integrate import quad

# The final equation for the initial sum of coordinates S(0) as a function of tau
# S(0; tau) = x(0; tau) + y(0; tau) + z(0; tau)
# We derived S(0; tau) = (3/4) * exp(2*tau**2) + 1/4
# The integral to be computed is the integral from 0 to infinity of 1/S(0; tau)

# Define the integrand function
def integrand(tau):
    """
    This is the function 1/S(0; tau) that we need to integrate.
    """
    s_0_tau = (3/4) * np.exp(2 * tau**2) + (1/4)
    return 1 / s_0_tau

# The final equation for the integral is:
# Integral from 0 to infinity of 4 / (3 * exp(2*tau**2) + 1) d(tau)
# We print the components of this equation as requested.
a = 4
b = 3
c = 2
d = 1
print(f"The integral to be evaluated is the integral from 0 to infinity of:")
print(f"{a} / ({b} * exp({c}*tau**2) + {d}) d(tau)")


# Perform the numerical integration from 0 to infinity
result, error = quad(lambda tau: 4 / (3 * np.exp(2 * tau**2) + 1), 0, np.inf)

print(f"\nThe numerical value of the integral is: {result}")
print(f"The analytical value of the integral is pi**2 / 8.")
print(f"The value of pi**2 / 8 is: {np.pi**2 / 8}")
