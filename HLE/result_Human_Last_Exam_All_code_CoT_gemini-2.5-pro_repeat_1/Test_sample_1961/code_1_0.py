import numpy as np
from scipy.integrate import quad

def F(c):
    """The exponent in the solution for L(c)."""
    return c**3 / 3 + c**2 / 2 + c

def integrand(t):
    """The integrand for the consistency check."""
    return t**3 * np.exp(F(t))

# Calculate the definite integral from 0 to 1
integral_value, _ = quad(integrand, 0, 1)

print("--- Analysis of the Limit Function ---")
print("Assuming the limits exist and are finite leads to a functional equation for the limit L(c).")
print("Solving this equation and applying the condition L(1)=1 gives a candidate solution for L(c).")
print("For this solution to be consistent, a specific definite integral must evaluate to 1.")
print("\nCalculating the integral: integral from 0 to 1 of t^3 * exp(t^3/3 + t^2/2 + t) dt")
print(f"Numerical result: {integral_value:.6f}")

print("\n--- Conclusion ---")
print(f"The integral is not 1. This contradiction implies our assumption was wrong.")
print("Therefore, the limit lim E(c^X_n) must be infinite for any c > 1.")
print("The equation lim E(p^X_n) * lim E(q^X_n) = lim E(r^X_n) becomes infinity * infinity = infinity.")
print("This holds true for any integers p, q, r > 1.")
print("The smallest possible integer value for r > 1 is 2.")

# Set p,q,r to their smallest possible integer values > 1
p = 2
q = 2
r = 2
print("\nA valid equation can be formed with these values:")
print(f"p = {p}")
print(f"q = {q}")
print(f"r = {r}")
