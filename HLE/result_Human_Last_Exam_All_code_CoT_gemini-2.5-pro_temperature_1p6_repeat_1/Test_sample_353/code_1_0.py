import numpy as np
import sympy as sp

def solve_bdf4_angle():
    """
    Finds the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Define the polynomial for x = cos(theta)
    # 18*x**3 - 9*x**2 - 11*x + 2 = 0
    coeffs = [18, -9, -11, 2]
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # The relevant root for the maximum angle is the largest real root in [-1, 1]
    # which is approximately 0.86438
    x0 = max(r.real for r in roots if np.isreal(r) and -1 <= r.real <= 1)
    
    # Expression for the real part of q(x)
    Re_q = 2*x0**4 - (16/3)*x0**3 + 4*x0**2 - 2/3
    
    # Expression for the imaginary part of q(x)
    # We take the positive square root for the angle in the first quadrant
    sin_theta = np.sqrt(1 - x0**2)
    Im_q_div_sin = -2*x0**3 + (16/3)*x0**2 - 5*x0 + 8/3
    Im_q = sin_theta * Im_q_div_sin
    
    # The angle alpha is arctan(Im(q)/Re(q))
    tan_alpha = Im_q / Re_q
    alpha_rad = np.arctan(tan_alpha)
    
    # --- Printing the exact value expression ---
    x_sym = sp.Symbol('x0')
    re_q_sym_str = f"2*({x_sym})**4 - (16/3)*({x_sym})**3 + 4*({x_sym})**2 - 2/3"
    im_q_sym_str = f"sqrt(1 - ({x_sym})**2) * (-2*({x_sym})**3 + (16/3)*({x_sym})**2 - 5*({x_sym}) + 8/3)"

    print("The A(alpha)-stability angle for BDF4 is given by alpha = arctan(Im(q)/Re(q)),")
    print("where q is evaluated at a specific theta_0, whose cosine x0 is a root of the polynomial:")
    print("18*x0**3 - 9*x0**2 - 11*x0 + 2 = 0")
    print(f"The relevant root is x0 = {x0:.10f}\n")
    
    print("The components of q(x0) are:")
    print(f"Re(q(x0)) = 2*({x0:.4f})**4 - (16/3)*({x0:.4f})**3 + 4*({x0:.4f})**2 - 2/3 = {Re_q:.10f}")
    print(f"Im(q(x0)) = sqrt(1 - ({x0:.4f})**2) * (-2*({x0:.4f})**3 + (16/3)*({x0:.4f})**2 - 5*({x0:.4f}) + 8/3) = {Im_q:.10f}\n")

    print("The final equation for the angle alpha (in radians) is:")
    print(f"alpha = arctan({Im_q:.10f} / {Re_q:.10f})")
    print(f"alpha = arctan({tan_alpha:.10f})")
    print(f"alpha = {alpha_rad:.10f}")

solve_bdf4_angle()