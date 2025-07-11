import numpy as np
import cmath

def solve_bdf4_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    
    # 1. Solve the cubic equation for c = cos(theta)
    # The equation is 50c^3 - 136c^2 + 123c - 40 = 0
    coeffs = [50, -136, 123, -40]
    roots = np.roots(coeffs)
    
    # Find the real root c0 in the range [-1, 1]
    c0 = 0
    for r in roots:
        if np.isreal(r) and -1 <= np.real(r) <= 1:
            c0 = np.real(r)
            break
            
    s0 = np.sqrt(1 - c0**2)
    
    # 2. Calculate the coordinates of the tangent point z(theta0)
    # x0 = 2*c0^4 - (16/3)*c0^3 + 4*c0^2 - 2/3
    # y0 = s0 * (-2*c0^3 + (16/3)*c0^2 - 5*c0 + 8/3)
    
    x0_num = 2*c0**4 - (16/3)*c0**3 + 4*c0**2 - (2/3)
    y0_num = s0 * (-2*c0**3 + (16/3)*c0**2 - 5*c0 + (8/3))
    
    # The stability angle alpha is arctan(-y0/x0)
    # Note: BDF stability region is in the left-half plane, so x0 should be negative.
    # The tangent point z0 is typically in the 2nd or 3rd quadrant.
    # We choose the tangent in the upper half-plane (y0>0), z0 in the 2nd quadrant (x0<0).
    # Then -z0 is in the 4th quadrant. The angle alpha is the angle of -z0 with the negative real axis.
    # alpha = arg(-z0) if -z0 is in the first quadrant, but here it's symmetric.
    # alpha = arctan(-y0 / x0) will give the angle in radians.

    # 3. Print the results in the required format.
    # We print the components to form the final expression for alpha.
    print(f"The cosine of the angle theta_0 is c0 = {c0}")
    print(f"The coordinates of the tangent point are (x0, y0) = ({x0_num}, {y0_num})")
    
    # We construct the final expression for alpha.
    # Since the problem asks for the *exact value* in terms of arctan(), we print the structure of the answer.
    # We output each number in the equation, as requested.
    
    y_val = -y0_num
    x_val = x0_num
    
    # Ensure a positive angle representation for the final answer.
    if y_val < 0:
        y_val = -y_val
        x_val = -x_val

    print("\nThe exact value of the angle alpha is given by the expression:")
    print(f"alpha = arctan({y_val} / {-x_val})")
    
    # We calculate the final numerical value for verification
    alpha_rad = np.arctan(-y0_num / x0_num)
    alpha_deg = np.rad2deg(alpha_rad)
    print(f"\nThis evaluates to approximately {alpha_rad} radians, or {alpha_deg} degrees.")

solve_bdf4_angle()