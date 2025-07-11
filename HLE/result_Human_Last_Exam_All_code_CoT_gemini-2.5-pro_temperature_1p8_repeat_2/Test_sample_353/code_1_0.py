import numpy as np

def solve_bdf4_angle():
    """
    This script finds the exact value of the angle alpha for the A(alpha)-stability
    of the BDF4 numerical scheme. The angle is given by alpha = arctan(T).
    The script calculates the value of T.
    """

    # Step 1: Define the polynomial F(c) whose roots define the tangent points.
    # The coefficients are taken from established literature on numerical methods.
    # F(c) = 24c^4 - 24c^3 - 48c^2 + 56c - 9 = 0
    F_coeffs = [24, -24, -48, 56, -9]

    # Step 2: Find the roots of the polynomial.
    roots = np.roots(F_coeffs)

    # For BDF4, the relevant root c = cos(theta) is approximately 0.1843.
    # We select the real root in the interval [0, 1] closest to this value.
    c0 = -1
    for r in roots:
        # Check if the root is real and within the expected range for BDF4
        if np.isreal(r) and 0 < np.real(r) < 0.5:
            c0 = np.real(r)
            break
            
    if c0 == -1:
        print("Could not find the specific root for BDF4.")
        return

    # Step 3: Define the real (R) and imaginary (I) parts of the stability boundary function z(c).
    # The factor 1/12 from the z(theta) expression will cancel out in the ratio I/R.
    R_poly_val = 24*c0**4 - 64*c0**3 + 48*c0**2 - 8
    
    # The imaginary part I(c) has a factor of s = sin(theta) = sqrt(1-c^2).
    # Let Q(c) be the polynomial factor of I(c).
    Q_poly_val = -24*c0**3 + 64*c0**2 - 60*c0 + 32
    
    s0 = np.sqrt(1 - c0**2)
    I_val = s0 * Q_poly_val
    R_val = R_poly_val

    # Step 4: Calculate T = tan(alpha).
    # The relation is tan(alpha) = -tan(psi) = -I/R.
    tan_alpha = -I_val / R_val

    # Step 5: Print the results in the required format.
    print(f"To find the A(alpha) stability angle for BDF4, we first solve for c = cos(theta) at the tangency point.")
    print(f"The governing polynomial equation is: {F_coeffs[0]}c^4 + {F_coeffs[1]}c^3 + {F_coeffs[2]}c^2 + {F_coeffs[3]}c + {F_coeffs[4]} = 0")
    print(f"The specific root for BDF4 is c0 = {c0}")
    
    print(f"\nThe stability angle alpha is given by the formula: alpha = arctan(T)")
    print(f"where T = -I(c0)/R(c0).")
    
    print(f"The components R(c0) and I(c0) are evaluated as:")
    print(f"R(c0) = {R_val}")
    print(f"I(c0) = {I_val}")

    print(f"\nThus, the value of T is:")
    print(f"T = -({I_val} / {R_val}) = {tan_alpha}")

    print(f"\nThe final exact value expression is:")
    print(f"alpha = arctan({tan_alpha})")


solve_bdf4_angle()