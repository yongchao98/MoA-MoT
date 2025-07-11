import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.
    """
    # 1. Define problem parameters
    # Let the chord length c be our reference unit.
    c = 1.0
    # Separation between trailing edge of aerofoil 1 and leading edge of aerofoil 2
    s = 0.5 * c
    # Ride height above the ground
    h = 0.5 * c

    print("--- Problem Setup ---")
    print(f"Chord c = {c}")
    print(f"Separation s = {s}")
    print(f"Ride Height h = {h}\n")

    # 2. Define the geometry of vortices and control points
    # Aerofoil 1 (front)
    # Vortex is at c/4, control point is at 3c/4 from the leading edge.
    P1_x, P1_y = c / 4.0, h
    C1_x, C1_y = 3.0 * c / 4.0, h

    # Aerofoil 2 (rear)
    # Its leading edge is at x = c + s
    P2_x, P2_y = (c + s) + c / 4.0, h
    C2_x, C2_y = (c + s) + 3.0 * c / 4.0, h

    # Mirror image vortices (circulation is opposite)
    P1_img_x, P1_img_y = P1_x, -h
    P2_img_x, P2_img_y = P2_x, -h

    # 3. Define the induced velocity function based on Biot-Savart law
    # The upwash `w` at a control point C from a vortex Γ at P is:
    # w = -Γ * (xC - xP) / (2 * pi * r^2)
    # We can define an influence coefficient C = (xC - xP) / r^2
    # The 1/(2*pi) factor cancels out in the final equation, so we can omit it.
    def get_influence_coefficient(vortex_pos, control_pos):
        """Calculates the influence coefficient C = (x_diff) / r^2."""
        x_diff = control_pos[0] - vortex_pos[0]
        y_diff = control_pos[1] - vortex_pos[1]
        r_sq = x_diff**2 + y_diff**2
        if r_sq == 0:
            return float('inf')
        return x_diff / r_sq

    # 4. Calculate influence coefficients for the equation
    # Equation form: w_total_1 = w_total_2
    # -(G2*C_21 - G1*C_1i1 - G2*C_2i1) = -(G1*C_12 - G1*C_1i2 - G2*C_2i2)
    # where C_ab is influence of vortex a on control point b

    # Influences on Control Point 1 (C1)
    C_21 = get_influence_coefficient((P2_x, P2_y), (C1_x, C1_y))
    C_1i1 = get_influence_coefficient((P1_img_x, P1_img_y), (C1_x, C1_y))
    C_2i1 = get_influence_coefficient((P2_img_x, P2_img_y), (C1_x, C1_y))

    # Influences on Control Point 2 (C2)
    C_12 = get_influence_coefficient((P1_x, P1_y), (C2_x, C2_y))
    C_1i2 = get_influence_coefficient((P1_img_x, P1_img_y), (C2_x, C2_y))
    C_2i2 = get_influence_coefficient((P2_img_x, P2_img_y), (C2_x, C2_y))
    
    print("--- Calculation ---")
    print("The flow tangency condition gives the equation:")
    print("Induced upwash at Control Point 1 = Induced upwash at Control Point 2")
    print("This leads to a linear equation for the circulations Gamma1 and Gamma2:\n")
    print("A * Gamma1 = B * Gamma2\n")

    # 5. Solve for the ratio Gamma1/Gamma2
    # Rearranging the equation:
    # Gamma1 * (-C_1i1 - C_12 + C_1i2) = Gamma2 * (-C_2i2 - C_21 + C_2i1)
    # A = (-C_1i1 - C_12 + C_1i2)
    # B = (-C_2i2 - C_21 + C_2i1)
    # Gamma1 / Gamma2 = B / A
    
    coeff_A = -C_1i1 - C_12 + C_1i2
    coeff_B = -C_2i2 - C_21 + C_2i1

    print("Where the coefficients are:")
    print(f"A = -({C_1i1:.3f}) - ({C_12:.3f}) + ({C_1i2:.3f}) = {coeff_A:.3f}")
    print(f"B = -({C_2i2:.3f}) - ({C_21:.3f}) + ({C_2i1:.3f}) = {coeff_B:.3f}\n")
    print(f"The equation is: {coeff_A:.3f} * Gamma1 = {coeff_B:.3f} * Gamma2")

    if coeff_A != 0:
        ratio = coeff_B / coeff_A
    else:
        ratio = float('inf')

    # 6. Final Result
    print("\nThe lift ratio L1/L2 is equal to the circulation ratio Gamma1/Gamma2.")
    print(f"L1/L2 = {coeff_B:.3f} / {coeff_A:.3f} = {ratio:.3f}")

calculate_lift_ratio()
<<< -0.2 >>>