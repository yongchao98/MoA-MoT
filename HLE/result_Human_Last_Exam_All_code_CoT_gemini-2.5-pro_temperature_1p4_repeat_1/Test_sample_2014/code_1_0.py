import numpy as np

def solve_aerofoil_tandem_ground_effect():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.
    """
    print("This script calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.")
    print("The problem is solved using the mirror image method and thin aerofoil theory.")
    print("\n--- Step 1: Define Model Geometry ---")

    # We can use dimensionless units by setting the chord c = 1.
    c = 1.0
    # Separation between aerofoils (Trailing Edge of aerofoil 1 to Leading Edge of aerofoil 2)
    s = 0.5 * c
    # Ride height above the ground
    h = 0.5 * c

    print(f"Chord length c = {c:.2f}")
    print(f"Separation s = 1/2 * c = {s:.2f}")
    print(f"Ride height h = 1/2 * c = {h:.2f}\n")

    # Using thin aerofoil theory, we place the lifting vortex at the quarter-chord (c/4)
    # and the control point (where flow tangency is enforced) at the three-quarter chord (3c/4).

    # Aerofoil 1 (leading)
    le1 = 0.0 # Leading edge of aerofoil 1 at x=0
    x_v1 = le1 + c / 4  # Position of vortex Gamma_1
    y_v1 = h
    x_p1 = le1 + 3 * c / 4 # Position of control point P1
    y_p1 = h

    # Aerofoil 2 (trailing)
    le2 = le1 + c + s # Leading edge of aerofoil 2
    x_v2 = le2 + c / 4  # Position of vortex Gamma_2
    y_v2 = h
    x_p2 = le2 + 3 * c / 4 # Position of control point P2
    y_p2 = h

    # Image vortices (to simulate the ground effect)
    # Vortices have opposite circulation (-Gamma) to satisfy the no-flow boundary condition at the ground.
    x_v1_img, y_v1_img = x_v1, -h
    x_v2_img, y_v2_img = x_v2, -h

    print("Vortex and Control Point Locations (c=1):")
    print(f"  Aerofoil 1: Vortex Γ1 at ({x_v1:.2f}, {y_v1:.2f}), Control Point P1 at ({x_p1:.2f}, {y_p1:.2f})")
    print(f"  Aerofoil 2: Vortex Γ2 at ({x_v2:.2f}, {y_v2:.2f}), Control Point P2 at ({x_p2:.2f}, {y_p2:.2f})")
    print(f"  Image 1:   Vortex -Γ1 at ({x_v1_img:.2f}, {y_v1_img:.2f})")
    print(f"  Image 2:   Vortex -Γ2 at ({x_v2_img:.2f}, {y_v2_img:.2f})")


    print("\n--- Step 2: Formulate System of Equations ---")
    # The vertical velocity (downwash) 'w' induced by a vortex 'Γ' at (vx, vy) on a point (px, py) is:
    # w = -Γ / (2*pi) * (px - vx) / r^2, where r^2 = (px - vx)^2 + (py - vy)^2.
    # The flow tangency equation is: Γ_i = A - π*c * w_i, where A = π*c*U*α is the isolated circulation.
    # w_i is the total downwash from all *other* vortices.

    def get_w_div_gamma(p_x, p_y, v_x, v_y):
        """Calculates the downwash per unit circulation (w/Γ) from a single vortex."""
        x_diff = p_x - v_x
        y_diff = p_y - v_y
        r_sq = x_diff**2 + y_diff**2
        if r_sq < 1e-9: return 0.0 # A vortex does not induce velocity at its own core
        return -1 / (2 * np.pi) * x_diff / r_sq

    # Downwash at P1 (w1) is induced by Γ2, -Γ1_img, and -Γ2_img.
    # w1 = (w/Γ from Γ2)*Γ2 + (w/Γ from -Γ1)*(-Γ1) + (w/Γ from -Γ2)*(-Γ2)
    w1_coeff_G1 = -get_w_div_gamma(x_p1, y_p1, x_v1_img, y_v1_img)
    w1_coeff_G2 = get_w_div_gamma(x_p1, y_p1, x_v2, y_v2) - get_w_div_gamma(x_p1, y_p1, x_v2_img, y_v2_img)

    # Downwash at P2 (w2) is induced by Γ1, -Γ1_img, and -Γ2_img.
    w2_coeff_G1 = get_w_div_gamma(x_p2, y_p2, x_v1, y_v1) - get_w_div_gamma(x_p2, y_p2, x_v1_img, y_v1_img)
    w2_coeff_G2 = -get_w_div_gamma(x_p2, y_p2, x_v2_img, y_v2_img)

    # Substitute w1 and w2 into the tangency equations to get a linear system for Γ1 and Γ2:
    # [ M11  M12 ] [Γ1] = [A]
    # [ M21  M22 ] [Γ2]   [A]
    M11 = 1 + np.pi * c * w1_coeff_G1
    M12 = np.pi * c * w1_coeff_G2
    M21 = np.pi * c * w2_coeff_G1
    M22 = 1 + np.pi * c * w2_coeff_G2

    print("\nThe resulting system of linear equations for circulations Γ1 and Γ2 is:")
    print(f"({M11:.3f}) * Γ1 + ({M12:.3f}) * Γ2 = A")
    print(f"({M21:.3f}) * Γ1 + ({M22:.3f}) * Γ2 = A")

    print("\n--- Step 3: Solve for the Lift Ratio ---")
    # Equating the two expressions for 'A':
    # M11*Γ1 + M12*Γ2 = M21*Γ1 + M22*Γ2
    # (M11 - M21)*Γ1 = (M22 - M12)*Γ2
    # Γ1/Γ2 = (M22 - M12) / (M11 - M21)

    numerator = M22 - M12
    denominator = M11 - M21
    lift_ratio = numerator / denominator

    print("By equating the two equations, we can solve for the ratio Γ1/Γ2.")
    print("The final equation is:")
    print(f"({denominator:.3f}) * Γ1 = ({numerator:.3f}) * Γ2")
    
    print("\nSince Lift L is proportional to Γ, the lift ratio L1/L2 equals Γ1/Γ2.")
    print(f"L1 / L2 = {numerator:.3f} / {denominator:.3f} = {lift_ratio:.3f}")
    
    # Store the final numerical answer in a specific format for parsing
    print(f"\n<<<ans>>>\n{lift_ratio}\n<</ans>>>")

solve_aerofoil_tandem_ground_effect()