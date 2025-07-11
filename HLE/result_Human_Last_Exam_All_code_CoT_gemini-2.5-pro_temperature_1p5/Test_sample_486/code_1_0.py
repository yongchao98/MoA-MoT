import numpy as np
from scipy.integrate import quad

# Step 1 & 2: Define the one-dimensional solution and its derivative.
# The PDE is Delta u = u^3 - u.
# A one-dimensional solution u(z) = tanh(z/sqrt(2)) solves this PDE and satisfies |u|<1.
#
# Step 3: Calculate |nabla u|^2 for this solution.
# For u(z), nabla u = (0, 0, u'(z)), so |nabla u|^2 = (u'(z))^2.
# The derivative u'(z) = d/dz [tanh(z/sqrt(2))] = (1/sqrt(2)) * sech^2(z/sqrt(2)).
# So, |nabla u|^2 = (1/2) * sech^4(z/sqrt(2)).

def integrand_sq_grad_u(z):
    """
    This function computes |nabla u|^2 = (u'(z))^2 for u(z) = tanh(z/sqrt(2)).
    Note that sech(x) = 1/cosh(x).
    """
    val = z / np.sqrt(2)
    sech_val = 1 / np.cosh(val)
    return 0.5 * (sech_val**4)

# Step 4: Analyze the integral over the ball B_R.
# The integral is I(R) = integral over B_R of |nabla u|^2 dV.
# We can evaluate this by slicing the ball into disks perpendicular to the z-axis.
# For a given z, the disk has radius sqrt(R^2 - z^2) and area pi*(R^2 - z^2).
# I(R) = integral from -R to R of [ pi*(R^2 - z^2) * |nabla u|^2 ] dz
# For large R, the term |nabla u|^2 (which is sech^4) decays very quickly,
# so the limits can be extended to +/- infinity.
# I(R) approx pi*R^2 * integral from -inf to inf of |nabla u|^2 dz
# Let C1 = integral from -inf to inf of |nabla u|^2 dz.
# We calculate C1 both numerically and analytically.

# Numerical calculation
C1_numerical, error = quad(integrand_sq_grad_u, -np.inf, np.inf)

# Analytical calculation
# The integral of 0.5 * sech^4(z/sqrt(2)) from -inf to inf is (2*sqrt(2))/3.
C1_analytical = 2 * np.sqrt(2) / 3

print("For the 1D solution u(z) = tanh(z/sqrt(2)), the integral of |nabla u|^2 over the ball B_R behaves as:")
print(f"Integral(|nabla u|^2, B_R) ≈ π * C1 * R^2 for large R.")
print(f"The constant C1 = ∫|∇u|² dz (from -∞ to ∞).")
print(f"Numerically calculated C1 = {C1_numerical:.8f}")
print(f"Analytically derived C1 = 2√2/3 ≈ {C1_analytical:.8f}")
print("-" * 30)

# Step 5: Determine the value of 'a'.
# We are given the condition: liminf_{R->∞} R^{-a} * Integral(|nabla u|^2, B_R) > 0
# Substituting the behavior of the integral:
# liminf_{R->∞} R^{-a} * (π * C1 * R^2) > 0
# This simplifies to: (π * C1) * liminf_{R->∞} R^(2-a) > 0
# Since π and C1 are positive constants, for the limit to be a positive constant,
# the exponent of R must be zero.

print("To satisfy the condition, we need the limit of R^(2-a) to be a positive constant.")
print("This leads to the final equation for 'a':")
final_equation_lhs_coeff = 2
final_equation_var = "a"
final_equation_rhs = 0
print(f"{final_equation_lhs_coeff} - {final_equation_var} = {final_equation_rhs}")

a = final_equation_lhs_coeff - final_equation_rhs
print(f"Solving for 'a', we get a = {a}.")
print("\nStep 6: General arguments (e.g., the monotonicity formula) show that this is the fastest possible growth rate for the integral.")
print("Therefore, the largest possible value for 'a' is 2.")
