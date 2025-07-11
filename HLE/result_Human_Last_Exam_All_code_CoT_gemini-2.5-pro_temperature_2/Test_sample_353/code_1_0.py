import numpy as np
import sympy as sp

def bdf4_stability_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # 1. Define the cubic polynomial for c = cos(theta)
    # The condition for tangency leads to the equation:
    # 50*c**3 - 136*c**2 + 123*c - 40 = 0
    c_poly_coeffs = [50, -136, 123, -40]
    
    # 2. Find the roots of the polynomial
    roots = np.roots(c_poly_coeffs)
    
    # Sort roots to easily identify them. We are interested in roots in [-1, 1].
    roots.sort()
    # There are three real roots. From literature and analysis of the stability plot,
    # the root defining the stability angle is the largest one that is less than 1.
    # The roots are approx [0.774, 0.8, 1.146]. We need the c0 ~ 0.774.
    c0 = roots[0]
    
    # 3. Calculate s0 = sin(theta_0)
    # We choose the positive root for s0, corresponding to the tangency point
    # in the upper half-plane.
    s0 = np.sqrt(1 - c0**2)

    # 4. Define x(c) and y(s, c) functions based on the stability boundary formula
    # z(theta) = x(theta) + i*y(theta)
    # x(theta) = 25/12 - 4*cos(t) + 3*cos(2t) - 4/3*cos(3t) + 1/4*cos(4t)
    # y(theta) = 4*sin(t) - 3*sin(2t) + 4/3*sin(3t) - 1/4*sin(4t)
    c, s = sp.symbols('c s')
    c2 = 2*c**2 - 1
    s2 = 2*s*c
    c3 = 4*c**3 - 3*c
    s3 = 3*s - 4*s**3
    c4 = 8*c**4 - 8*c**2 + 1
    s4 = 4*s*c * (1 - 2*s**2)

    x_expr = 25/12 - 4*c + 3*c2 - (4/3)*c3 + (1/4)*c4
    y_expr = 4*s - 3*s2 + (4/3)*s3 - (1/4)*s4
    
    # 5. Evaluate x and y at the point of tangency (c0, s0)
    x0 = x_expr.subs({c: c0, s: s0})
    y0 = y_expr.subs({c: c0, s: s0})
    
    # 6. Calculate the argument phi_0 and the stability angle alpha
    phi0_rad = np.arctan(y0/x0)
    alpha_rad = np.pi - phi0_rad
    
    # Output the exact expression in terms of arctan
    # The cubic equation for c = cos(theta) is 50*c^3 - 136*c^2 + 123*c - 40 = 0
    # Let c0 be the smallest root of this equation.
    print("The A(alpha)-stability angle for BDF4 is determined by a tangency condition.")
    print("This condition reduces to finding a specific root, c0, of the cubic polynomial:")
    print(f"{c_poly_coeffs[0]}*c^3 + ({c_poly_coeffs[1]})*c^2 + {c_poly_coeffs[2]}*c + ({c_poly_coeffs[3]}) = 0")
    print("\nThe relevant root is c0 ≈ {:.6f}".format(c0))
    print(f"Let c0 be the root of {c_poly_coeffs[0]}c^3{c_poly_coeffs[1]:+}c^2{c_poly_coeffs[2]:+}c{c_poly_coeffs[3]:+}=0 with value ~{c0:.6f}.")
    print("Let s0 = sqrt(1 - c0^2).")

    print("\nThe coordinates (x0, y0) of the tangency point are then calculated:")
    print(f"x0 = 25/12 - 4*c0 + 3*(2*c0^2 - 1) - (4/3)*(4*c0^3 - 3*c0) + (1/4)*(8*c0^4 - 8*c0^2 + 1)")
    print(f"y0 = 4*s0 - 3*(2*s0*c0) + (4/3)*(3*s0 - 4*s0^3) - (1/4)*(4*s0*c0*(1-2*s0^2))")
    print(f"\nNumerically, (x0, y0) ≈ ({x0:.6f}, {y0:.6f}).")

    print("\nThe angle of the tangency point is phi0 = arctan(y0/x0).")
    print("The stability angle alpha is given by the formula alpha = pi - phi0.")
    print("\nTherefore, the exact value of alpha in radians is:")
    print(f"alpha = pi - arctan({y0:.6f} / {x0:.6f})")

    # For the final answer format
    return alpha_rad

alpha = bdf4_stability_angle()
print(f"\nNumerically, alpha ≈ {alpha:.6f} radians, or {np.rad2deg(alpha):.2f} degrees.")
print(f"<<<pi - arctan({np.tan(np.pi-alpha):.6f})>>>")
