import numpy as np
from scipy.integrate import quad

# This script provides numerical evidence to determine the value of 'a'.
# The problem is to find the largest 'a' such that the following limit is positive:
# liminf_{R->inf} R^{-a} * integral_{B_R} |nabla u|^2 > 0
# for solutions to Delta u = u^3 - u in R^3 with |u|<1.
#
# Our analysis suggests that the maximal growth of the Dirichlet integral (integral_{B_R} |nabla u|^2)
# occurs for one-dimensional solutions, where the integral grows proportionally to R^2. This implies a=2.
#
# A well-known 1D solution is u(x1, x2, x3) = tanh(x3 / sqrt(2)).
# For this solution, the squared gradient is |nabla u|^2 = (1/2) * sech^4(x3 / sqrt(2)).
#
# To compute the integral over a ball B_R, we can integrate with respect to x3 from -R to R.
# At each height x3, the cross-section of the ball is a disk of radius sqrt(R^2 - x3^2) and
# area pi*(R^2 - x3^2).
#
# The full integral is I(R) = integral from -R to R of [ pi*(R^2 - z^2) * (1/2)*sech^4(z/sqrt(2)) ] dz.
# We will numerically compute this integral for several large values of R and show that the ratio I(R)/R^2
# converges to a constant. This provides strong numerical evidence that a=2.

def integrand(z, R):
    """
    This is the function to be integrated with respect to z from -R to R.
    It represents the contribution to the total integral from the disk at height z.
    """
    # sech(x) is the hyperbolic secant, equal to 1/cosh(x)
    sech_val = 1.0 / np.cosh(z / np.sqrt(2))
    return np.pi * (R**2 - z**2) * 0.5 * (sech_val**4)

# Calculate the integral for a sequence of increasing radii R and compute the ratio I(R)/R^2.
print("Calculating the integral I(R) = integral_{B_R} |nabla u|^2 for increasing values of R.")
print("If 'a' were 2, the ratio I(R) / R^2 should converge to a non-zero constant.")
print()
# The final "equation" is I(R) ~ C * R^a. We test a=2 and find C.
print(f"{'R':<10} | {'I(R)':<25} | {'C = I(R)/R^2':<25}")
print("-" * 65)

# A list of radii to test. They should be large enough for the asymptotic behavior to appear.
R_values = [10, 20, 50, 100, 200, 500]

for R in R_values:
    # Perform the numerical integration from -R to R using SciPy's quad function.
    # The 'args' parameter passes the current value of R to the integrand.
    integral_value, error_estimate = quad(integrand, -R, R, args=(R,))

    # Calculate the ratio to test the R^2 growth hypothesis.
    ratio = integral_value / (R**2)

    # Print the results for this R in a formatted table.
    # We output each relevant number: R, the integral value, and the computed constant C.
    print(f"{R:<10.1f} | {integral_value:<25.10f} | {ratio:<25.10f}")

# The constant C that the ratio is converging to can also be computed analytically.
# As R -> infinity, I(R) ~ (pi/2) * R^2 * integral_{-inf}^{+inf} sech^4(z/sqrt(2)) dz.
# Let y = z/sqrt(2), so dz = sqrt(2) dy. The definite integral evaluates to (4/3).
# C_limit = (pi/2) * sqrt(2) * (4/3) = (2*sqrt(2)/3)*pi.
analytic_C = (2 * np.sqrt(2) / 3) * np.pi
print("-" * 65)
print(f"The theoretical limit for the ratio C is (2*sqrt(2)/3)*pi approx {analytic_C:.10f}")
print("\nAs the table clearly shows, the ratio I(R)/R^2 converges to this theoretical constant.")
print("This provides strong evidence that the Dirichlet integral grows quadratically with R.")
print("Therefore, the largest possible value for 'a' is 2.")
