import math

def solve_beam_force():
    """
    This function calculates the force F required to make the deflection at the
    end of the beam zero.
    """
    # Step 1: Calculate I_zz and I_ss
    # Given a = 12^(1/4), so a^4 = 12.
    a_fourth_power = 12.0

    # For a square b x b, the moment of inertia about a centroidal axis is b^4 / 12.
    # Main square (3a x 3a):
    I_main_c = (3**4 * a_fourth_power) / 12
    
    # Cutout squares (a x a):
    I_cutout_c = a_fourth_power / 12
    Area_cutout = a_fourth_power**0.5 * a_fourth_power**0.5 # which is a^2

    # Using the parallel axis theorem: I = I_c + A * d^2
    # I_zz = integral(s^2 dA). The distance 'd' is the s-coordinate of the cutout's center.
    # The centers are at s = a/2 and s = -a/2.
    # d_s_sq = (a/2)^2 = a^2 / 4 = sqrt(a^4) / 4
    d_s_sq = math.sqrt(a_fourth_power) / 4 # This is incorrect. d^2 = (a/2)^2 = a^2/4. a^2 = sqrt(a^4).
    a_sq = math.sqrt(a_fourth_power)
    d_s_sq = (a_sq / 2)**2

    # Let's use the derived formulas to avoid potential floating point issues with roots
    # I_zz = I_main_c - 2 * (I_cutout_c + Area_cutout * d_s^2)
    # I_zz = (81*a^4/12) - 2 * (a^4/12 + a^2 * (a/2)^2)
    # I_zz = (81*a^4/12) - 2 * (a^4/12 + a^4/4) = (81*a^4/12) - 2 * (4*a^4/12)
    # I_zz = (81*a^4/12) - (8*a^4/12) = 73*a^4/12
    I_zz = 73 * a_fourth_power / 12

    # I_ss = integral(z^2 dA). The distance 'd' is the z-coordinate of the cutout's center.
    # The centers are at z = -a and z = a. d_z_sq = a^2.
    # I_ss = I_main_c - 2 * (I_cutout_c + Area_cutout * d_z^2)
    # I_ss = (81*a^4/12) - 2 * (a^4/12 + a^2 * a^2)
    # I_ss = (81*a^4/12) - 2 * (a^4/12 + a^4) = (81*a^4/12) - 2 * (13*a^4/12)
    # I_ss = (81*a^4/12) - (26*a^4/12) = 55*a^4/12
    I_ss = 55 * a_fourth_power / 12

    # Step 2: Calculate L and q_0 from the given data
    L = 30 * I_zz / 73
    q_0 = 9 * I_ss / 55
    
    # Step 3: Use the derived formula to calculate F
    # The derivation y_total(3L/2) = 0 leads to:
    # 9*F*L^3 / (8*EI) = 13*q_0*L^4 / (240*EI)
    # F = (13 * q_0 * L * 8) / (240 * 9)
    # F = 13 * q_0 * L / 270
    F = (13 * q_0 * L) / 270

    # Step 4: Print the equation with numerical values and the final answer
    print("The final equation for F is derived by setting the total deflection to zero.")
    print("Deflection from F + Deflection from q(x) = 0")
    print(f"This results in the equation: F = (13 * q_0 * L) / 270")
    print("\nSubstituting the calculated values for q_0 and L:")
    print(f"F = (13 * {q_0:.2f} * {L:.2f}) / 270")
    
    print("\nThe calculated force F is:")
    print(f"{F:.2f}")

solve_beam_force()