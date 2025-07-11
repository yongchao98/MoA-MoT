import numpy as np

def solve_bdf4_angle():
    """
    Finds the exact value of the angle alpha for the A(alpha)-stability
    of the BDF4 scheme.
    """
    # Step 1: Define the polynomial for c = cos(theta)
    # The polynomial is 8c^3 + 24c^2 - 20c - 17 = 0
    p = np.poly1d([8, 24, -20, -17])
    
    # Step 2: Find the roots of the polynomial
    roots = p.roots
    
    # Filter for real roots in the valid range [-1, 1] for cosine
    c_roots = [r.real for r in roots if np.isreal(r) and -1 <= r.real <= 1]
    
    # The root that defines the angle alpha for BDF4 is the negative one.
    c0 = min(c_roots)
    
    # Step 3: Calculate the corresponding theta
    # We take the value in [0, pi]
    theta0 = np.arccos(c0)
    
    # Step 4: Evaluate z(theta0) on the stability boundary
    # z(theta) = (25 - 48*exp(-i*theta) + 36*exp(-i*2*theta) - 16*exp(-i*3*theta) + 3*exp(-i*4*theta))/12
    eta = np.exp(-1j * theta0)
    z0 = (25 - 48*eta + 36*eta**2 - 16*eta**3 + 3*eta**4) / 12
    x0 = z0.real
    y0 = z0.imag
    
    # The A(alpha)-stability angle alpha is given by pi - |arg(z0)|
    # Since z0 is in the second quadrant (x0 < 0, y0 > 0), arg(z0) = atan(y0/x0) + pi
    # So alpha = pi - (atan(y0/x0) + pi) = -atan(y0/x0)
    # Since alpha must be positive, alpha = atan(-y0/x0) = atan(y0/|x0|)
    
    # For printing, let's use the fraction.
    num = y0
    den = -x0
    
    print("The A(alpha)-stability angle for BDF4 is alpha, where:")
    # We express alpha in terms of arctan of the ratio of y0 and |x0|.
    # These values depend on the specific root of the cubic polynomial.
    print(f"alpha = arctan({num} / {den})")
    print(f"alpha = arctan({num/den})")
    alpha_rad = np.arctan(num/den)
    print(f"alpha ≈ {alpha_rad} radians")
    print(f"alpha ≈ {np.rad2deg(alpha_rad)} degrees")
    
    # To express the answer "exactly" we can use the polynomials for x0 and y0
    # in terms of c0=cos(theta0) and s0=sin(theta0).
    # From derivation:
    # 12*x0 = 24*c0**4 - 64*c0**3 + 48*c0**2 - 8
    # 12*y0 = s0 * (-24*c0**3 + 64*c0**2 - 60*c0 + 32)
    # where c0 is the specific root of 8c^3 + 24c^2 - 20c - 17 = 0.
    
    print("\nThe exact value is of the form:")
    print("alpha = arctan(Y/X) where X and Y are functions of c0, the root of 8c^3 + 24c^2 - 20c - 17 = 0.")
    print("Specifically, c0 ≈ {:.6f}".format(c0))
    print("alpha = arctan(sqrt(1-c0^2) * (-24*c0^3 + 64*c0^2 - 60*c0 + 32) / -(24*c0^4 - 64*c0^3 + 48*c0^2 - 8))")


solve_bdf4_angle()