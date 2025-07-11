import fractions

def solve_aerofoil_problem():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.
    """
    # Use fractions for exact arithmetic
    Fraction = fractions.Fraction

    # --- Step 1: Define the problem setup ---
    # Given: separation s = 1/2 * c and ride height h = 1/2 * c.
    # We can work with non-dimensional lengths by setting c=1.
    s = Fraction(1, 2)
    h = Fraction(1, 2)

    # Standard thin aerofoil model:
    # Vortex placed at quarter-chord (c/4)
    # Control point (where flow tangency is enforced) at three-quarter chord (3/4c)
    x_vortex_pos = Fraction(1, 4)
    x_control_point_pos = Fraction(3, 4)

    # --- Step 2: Define the locations of vortices and control points ---
    # Aerofoil 1 (leading)
    x_v1 = x_vortex_pos
    y_v1 = h
    x_cp1 = x_control_point_pos
    y_cp1 = h

    # Aerofoil 2 (trailing)
    x_v2 = s + x_vortex_pos
    y_v2 = h
    x_cp2 = s + x_control_point_pos
    y_cp2 = h

    # Image Aerofoil 1
    x_v1_img = x_v1
    y_v1_img = -h

    # Image Aerofoil 2
    x_v2_img = x_v2
    y_v2_img = -h

    # --- Step 3: Calculate the downwash coefficients ---
    # The vertical velocity (downwash) w induced at a point (xp, yp) by a vortex
    # of strength Gamma at (xv, yv) is: w = Gamma * (xp - xv) / (2 * pi * r^2)
    # We define a non-dimensional downwash coefficient k such that w = k * Gamma / (pi * c)
    # k = (xp - xv) / (2 * r^2)
    def get_downwash_coeff(xp, yp, xv, yv):
        """Calculates the non-dimensional downwash coefficient k."""
        x_rel = xp - xv
        y_rel = yp - yv
        if x_rel == 0 and y_rel == 0:
            return Fraction(0)
        r_sq = x_rel**2 + y_rel**2
        if r_sq == 0:
            return Fraction(0)
        return x_rel / (2 * r_sq)

    # Coefficients for downwash w1 at Control Point 1 (CP1)
    # w1 is induced by V2, V1_img, V2_img
    # w1/(pi*c) = k11*Gamma1 + k12*Gamma2
    k11 = -get_downwash_coeff(x_cp1, y_cp1, x_v1_img, y_v1_img) # From -Gamma1 at V1'
    k12 = get_downwash_coeff(x_cp1, y_cp1, x_v2, y_v2) - get_downwash_coeff(x_cp1, y_cp1, x_v2_img, y_v2_img)

    # Coefficients for downwash w2 at Control Point 2 (CP2)
    # w2 is induced by V1, V1_img, V2_img
    # w2/(pi*c) = k21*Gamma1 + k22*Gamma2
    k21 = get_downwash_coeff(x_cp2, y_cp2, x_v1, y_v1) - get_downwash_coeff(x_cp2, y_cp2, x_v1_img, y_v1_img)
    k22 = -get_downwash_coeff(x_cp2, y_cp2, x_v2_img, y_v2_img) # From -Gamma2 at V2'

    # --- Step 4: Set up and solve the system of linear equations ---
    # The flow tangency condition gives: Gamma = Gamma0 - pi * c * w
    # Eq1: G1 = G0 - (k11*G1 + k12*G2) => (1 + k11)*G1 + k12*G2 = G0
    # Eq2: G2 = G0 - (k21*G1 + k22*G2) => k21*G1 + (1 + k22)*G2 = G0

    # Solve for G1/G0 from Eq1 (since k12=0, it simplifies)
    G1_over_G0 = 1 / (1 + k11)

    # Substitute G1 into Eq2 and solve for G2/G0
    # k21*(G1_over_G0*G0) + (1 + k22)*G2 = G0
    # (1 + k22)*G2 = G0 * (1 - k21*G1_over_G0)
    G2_over_G0 = (1 - k21 * G1_over_G0) / (1 + k22)

    # --- Step 5: Calculate the final lift ratio ---
    # L1/L2 = Gamma1/Gamma2 = (Gamma1/Gamma0) / (Gamma2/Gamma0)
    lift_ratio = G1_over_G0 / G2_over_G0

    # --- Step 6: Print the results ---
    print("This script calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.")
    print(f"Configuration: Separation s = {float(s):.1f}c, Ride Height h = {float(h):.1f}c\n")
    print("--- Calculation Steps ---")
    print("1. The circulation of the first aerofoil relative to an isolated aerofoil (Gamma0) is:")
    print(f"   Gamma1/Gamma0 = 1 / (1 + ({k11})) = {G1_over_G0.numerator}/{G1_over_G0.denominator}")
    print("\n2. The circulation of the second aerofoil relative to an isolated aerofoil is:")
    print(f"   Gamma2/Gamma0 = (1 - ({k21})*({G1_over_G0})) / (1 + ({k22})) = {G2_over_G0.numerator}/{G2_over_G0.denominator}")
    print("\n3. The final lift ratio is the ratio of these two values:")
    print(f"   L1/L2 = (Gamma1/Gamma0) / (Gamma2/Gamma0)")
    print(f"   L1/L2 = ({G1_over_G0.numerator}/{G1_over_G0.denominator}) / ({G2_over_G0.numerator}/{G2_over_G0.denominator})")
    print(f"   L1/L2 = {lift_ratio.numerator}/{lift_ratio.denominator}")
    print(f"\n   The final result is L1/L2 ≈ {float(lift_ratio):.4f}")

if __name__ == '__main__':
    solve_aerofoil_problem()