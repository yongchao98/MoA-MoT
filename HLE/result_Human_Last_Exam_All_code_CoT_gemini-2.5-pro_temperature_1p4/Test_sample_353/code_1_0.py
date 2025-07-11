import numpy as np
import cmath

def bdf4_alpha_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Define and solve the quartic equation for c = cos(theta)
    # The equation for the tangency condition is:
    # 50*c^4 + 452*c^3 - 64*c^2 - 412*c + 63 = 0
    coeffs = [50, 452, -64, -412, 63]
    roots = np.roots(coeffs)

    max_alpha = 0
    tangent_val = 0
    
    # Step 2: Iterate through real roots c in [-1, 1]
    for c in roots:
        if np.isreal(c) and -1 <= c.real <= 1:
            c = c.real
            
            # We consider two possible angles theta for each c, one in [0, pi] and one in [pi, 2pi]
            # Since the stability region is symmetric with respect to the real axis,
            # we only need to consider theta in [0, pi], which means s >= 0.
            theta = np.arccos(c)
            s = np.sin(theta)

            # Step 3: Calculate the point z(theta) on the stability boundary
            # z(theta) = x + iy
            # c_k = cos(k*theta), s_k = sin(k*theta)
            c2 = np.cos(2*theta)
            s2 = np.sin(2*theta)
            c3 = np.cos(3*theta)
            s3 = np.sin(3*theta)
            c4 = np.cos(4*theta)
            s4 = np.sin(4*theta)
            
            x = 25/12 - 4*c + 3*c2 - (4/3)*c3 + (1/4)*c4
            y = 4*s - 3*s2 + (4/3)*s3 - (1/4)*s4
            
            # The stability region for BDF4 is the exterior of the curve z(theta)
            # and it extends into the left-half plane. The tangent from the origin
            # will be to a point in the LHP (x < 0).
            if x < 0:
                # Step 4: Calculate alpha
                # arg(z) is in (-pi, pi]. We use atan2 for correctness.
                arg_z = np.arctan2(y, x)
                # alpha must be positive, corresponds to angle in upper-half plane.
                if arg_z > np.pi/2:
                    alpha = np.pi - arg_z
                    if alpha > max_alpha:
                        max_alpha = alpha
                        # We want tan(alpha)
                        tangent_val = np.tan(alpha)

    # Print the equation in the required format
    alpha_deg = np.rad2deg(max_alpha)
    print(f"The equation for the stability angle alpha is alpha = arctan(T), where T is the tangent of alpha.")
    print(f"The numerical value for T = tan(alpha) is approximately {tangent_val:.8f}")
    print(f"The angle alpha is approximately {alpha_deg:.8f} degrees.")
    print("\nThe final answer expressed in the required format is:")
    print(f"alpha = arctan({tangent_val})")


if __name__ == '__main__':
    bdf4_alpha_angle()
