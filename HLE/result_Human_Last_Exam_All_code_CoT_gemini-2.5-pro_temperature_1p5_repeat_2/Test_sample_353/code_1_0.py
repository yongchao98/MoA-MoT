import numpy as np
from numpy.polynomial import Polynomial

# The A(alpha)-stability angle for BDF4 is determined by the point of tangency
# on its stability boundary. This leads to a specific value for cos(theta),
# which is a root of a cubic equation.

# The derived cubic polynomial for c = 2*cos(theta) is:
# 15*c^3 - 159*c^2 + 462*c - 404 = 0
coeffs = [-404, 462, -159, 15]
cubic_poly = Polynomial(coeffs)

# Find the roots of the cubic polynomial
roots = cubic_poly.roots()

# The physically meaningful solution for theta corresponds to a real root c
# such that c/2 is in [-1, 1].
valid_c = None
for r in roots:
    # Check if the root is real (within a small tolerance for numerical errors)
    if np.isclose(np.imag(r), 0):
        c = np.real(r)
        # c = 2*cos(theta), so c must be in [-2, 2]
        if -2 <= c <= 2:
            valid_c = c
            break

if valid_c is None:
    print("Could not find a valid root for c = 2*cos(theta).")
else:
    # From c, find cos(theta) and sin(theta) for the tangency point
    cos_theta = valid_c / 2.0
    # sin(theta) can be positive or negative; due to symmetry we can choose positive
    sin_theta = np.sqrt(1.0 - cos_theta**2)

    # The point on the unit circle is xi = cos(theta) + i*sin(theta)
    xi = cos_theta + 1j * sin_theta

    # Evaluate the stability function z(xi) = rho(xi) / xi^4 at this point.
    # It is numerically more stable to evaluate z(xi) using the form:
    # z(xi) = 25/12 - 4/xi + 3/xi^2 - 4/(3*xi^3) + 1/(4*xi^4)
    z_val = (25.0/12.0) - 4.0/xi + 3.0/xi**2 - (4.0/3.0)/xi**3 + (1.0/4.0)/xi**4

    # The stability angle alpha is given by tan(alpha) = -Im(z) / Re(z)
    tan_alpha = -np.imag(z_val) / np.real(z_val)
    alpha_rad = np.arctan(tan_alpha)

    # The final equation is alpha = arctan(T)
    print("The equation for the stability angle is alpha = arctan(T), where T is:")
    print(f"T = {tan_alpha}")
    print("\nThis is derived from the root of the cubic equation:")
    print(f"{coeffs[3]}c^3 + ({coeffs[2]})c^2 + {coeffs[1]}c + ({coeffs[0]}) = 0")
    print(f"The angle alpha is therefore:")
    print(f"alpha = arctan({tan_alpha:.10f}) = {alpha_rad:.10f} radians.")

<<<alpha = arctan(2.9852226245)>>>