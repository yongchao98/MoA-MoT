import math

def solve_aerofoil_tandem_ground_effect():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    formation in ground effect using the mirror image method.
    """

    # --- 1. Define parameters ---
    # We can set the chord c = 1.0, as the final ratio is non-dimensional.
    c = 1.0
    # Separation between trailing edge of aerofoil 1 and leading edge of aerofoil 2
    s = 0.5 * c
    # Ride height above the ground
    h = 0.5 * c

    print("--- Problem Setup ---")
    print(f"Chord length, c = {c}")
    print(f"Separation, s = {s}")
    print(f"Ride height, h = {h}\n")


    # --- 2. Define geometry (quarter-chord vortex, 3/4-chord control point model) ---
    # Aerofoil 1 (leading)
    x_vortex1 = c / 4.0
    z_vortex1 = h
    x_cp1 = 3.0 * c / 4.0
    z_cp1 = h

    # Aerofoil 2 (trailing)
    x_vortex2 = c + s + (c / 4.0)
    z_vortex2 = h
    x_cp2 = c + s + (3.0 * c / 4.0)
    z_cp2 = h

    # Mirror image vortices (for ground effect). Vortex strength is opposite (-Gamma).
    x_mirror_vortex1 = x_vortex1
    z_mirror_vortex1 = -h
    x_mirror_vortex2 = x_vortex2
    z_mirror_vortex2 = -h


    # --- 3. Define function to calculate induced velocity coefficient ---
    # The vertical velocity 'w' induced by a vortex Gamma is w = Gamma * C_w, where
    # C_w = -(1 / (2*pi)) * (delta_x / r^2)
    def get_w_coeff(x_v, z_v, x_p, z_p):
        """Calculates the coefficient for vertical induced velocity."""
        delta_x = x_p - x_v
        delta_z = z_p - z_v
        r_sq = delta_x**2 + delta_z**2
        if r_sq == 0:
            return 0  # Should not happen in this problem setup
        return -(1.0 / (2.0 * math.pi)) * (delta_x / r_sq)


    # --- 4. Set up the system of linear equations ---
    # The flow tangency condition gives two equations:
    # Eq1: C11*Gamma1 + C12*Gamma2 = K' (where K' is related to angle of attack)
    # Eq2: C21*Gamma1 + C22*Gamma2 = K'
    # We need to find the influence coefficients Cij.

    # Row 1: Induced velocity at Control Point 1 (CP1)
    # C11: Effect of Aerofoil 1's image vortex on CP1
    C11 = get_w_coeff(x_mirror_vortex1, z_mirror_vortex1, x_cp1, z_cp1) * (-1)
    # C12: Effect of Aerofoil 2's vortex and its image on CP1
    C12 = get_w_coeff(x_vortex2, z_vortex2, x_cp1, z_cp1) + \
          get_w_coeff(x_mirror_vortex2, z_mirror_vortex2, x_cp1, z_cp1) * (-1)

    # Row 2: Induced velocity at Control Point 2 (CP2)
    # C21: Effect of Aerofoil 1's vortex and its image on CP2
    C21 = get_w_coeff(x_vortex1, z_vortex1, x_cp2, z_cp2) + \
          get_w_coeff(x_mirror_vortex1, z_mirror_vortex1, x_cp2, z_cp2) * (-1)
    # C22: Effect of Aerofoil 2's image vortex on CP2
    C22 = get_w_coeff(x_mirror_vortex2, z_mirror_vortex2, x_cp2, z_cp2) * (-1)

    print("--- System of Equations ---")
    print("The system is of the form:")
    print("C11 * Γ1 + C12 * Γ2 = K")
    print("C21 * Γ1 + C22 * Γ2 = K")
    print(f"Calculated influence coefficients (Cij):")
    print(f"C11 = {C11:.4f}, C12 = {C12:.4f}")
    print(f"C21 = {C21:.4f}, C22 = {C22:.4f}\n")


    # --- 5. Solve the system for the ratio Gamma1 / Gamma2 ---
    # C11*Γ1 + C12*Γ2 = C21*Γ1 + C22*Γ2
    # (C11 - C21)*Γ1 = (C22 - C12)*Γ2
    # Γ1 / Γ2 = (C22 - C12) / (C11 - C21)

    numerator = C22 - C12
    denominator = C11 - C21
    
    # Check for division by zero
    if denominator == 0:
        if numerator == 0:
            ratio = float('nan') # Indeterminate
        else:
            ratio = float('inf') # Infinite
    else:
        ratio = numerator / denominator


    # --- 6. Print the final calculation and result ---
    print("--- Final Calculation ---")
    print("The lift ratio L1/L2 is equal to the vortex strength ratio Γ1/Γ2.")
    print("Γ1 / Γ2 = (C22 - C12) / (C11 - C21)")
    print(f"Γ1 / Γ2 = ({C22:.4f} - {C12:.4f}) / ({C11:.4f} - {C21:.4f})")
    print(f"Γ1 / Γ2 = ({numerator:.4f}) / ({denominator:.4f})")
    print(f"Γ1 / Γ2 = {ratio:.4f}\n")
    
    print("--- Result ---")
    print(f"The calculated lift ratio L1/L2 is: {ratio:.4f}")
    return ratio

# Execute the function
final_ratio = solve_aerofoil_tandem_ground_effect()
# The final answer in the required format
# <<< [final_ratio] >>>