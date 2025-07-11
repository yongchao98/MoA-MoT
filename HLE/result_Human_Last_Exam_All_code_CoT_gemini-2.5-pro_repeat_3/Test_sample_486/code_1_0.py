import numpy as np
from scipy.integrate import quad

# The problem is to find the largest 'a' such that for a solution u to Delta u = u^3 - u,
# the quantity liminf_{R->inf} R^{-a} * integral_{B_R} |nabla u|^2 > 0.

# This is equivalent to finding the maximum possible asymptotic growth rate of the integral.
# Theory suggests that for solutions representing stable minimal surfaces (planes),
# the growth rate of the integral is proportional to the area of the surface, which is O(R^2).
# This implies a <= 2.

# We verify that a=2 is achievable using the 1D solution u(z) = tanh(z/sqrt(2)).
# For this solution, |nabla u|^2 = (1/2) * sech^4(z/sqrt(2)).
# The integral over the ball B_R is approximately C * R^2 for large R, where
# C = (pi/2) * integral_{-inf}^{+inf} sech^4(z/sqrt(2)) dz.

# We calculate this constant C.

# Define the sech function
def sech(x):
    return 1 / np.cosh(x)

# Define the integrand from the asymptotic analysis
def integrand(z):
    return sech(z / np.sqrt(2))**4

# The integral of sech^4(x) from -inf to inf is 4/3.
# By substitution s = z/sqrt(2), the integral of sech^4(z/sqrt(2)) dz is (4/3)*sqrt(2).
analytic_integral_val = 4 * np.sqrt(2) / 3

# Perform numerical integration to verify
numeric_integral_val, error = quad(integrand, -np.inf, np.inf)

# The coefficient of R^2 is (pi/2) times this integral value
coefficient = (np.pi / 2) * numeric_integral_val

# The analysis shows the integral grows as R^a.
a = 2

print(f"For the 1D solution u(z) = tanh(z/sqrt(2)), the integral ∫|∇u|² over B_R has an asymptotic behavior.")
print(f"The integral is approximately C * R^a for large R.")
print(f"The power 'a' is determined by the dimensionality of the interface, which is a plane (2D) in a 3D space.")
print(f"Therefore, the power a = {a}")
print(f"The coefficient C = (π/2) * ∫sech⁴(z/√2)dz can be calculated numerically.")
print(f"Value of the integral part: {numeric_integral_val:.4f} (Analytic: {analytic_integral_val:.4f})")
print(f"The full coefficient C is (π/2) * {numeric_integral_val:.4f} ≈ {coefficient:.4f}")
print(f"Since a solution exists for which the integral grows as R^{a}, with a=2, and theory suggests no faster growth is possible, the largest possible value for 'a' is 2.")
