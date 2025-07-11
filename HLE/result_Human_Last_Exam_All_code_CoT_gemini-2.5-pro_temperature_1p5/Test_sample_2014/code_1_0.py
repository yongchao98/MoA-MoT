import math

def solve_tandem_aerofoil_ground_effect():
    """
    Calculates the lift ratio L1/L2 for two tandem aerofoils in ground effect.
    """
    # Step 1: Explain the model and define parameters
    print("Step 1: Define the geometry of the system.")
    print("Let the chord length c = 1 for simplicity.")
    c = 1.0
    s = 0.5 * c  # Separation between aerofoils
    h = 0.5 * c  # Ride height above the ground
    print(f"Chord c = {c}, Separation s = {s}, Ride Height h = {h}\n")

    print("Step 2: Define vortex coordinates based on the geometry.")
    # Real vortices (representing the aerofoils)
    vortex1_pos = (0, h)
    vortex2_pos = (s, h)
    # Image vortices (representing the ground effect)
    image1_pos = (0, -h)
    image2_pos = (s, -h)
    print(f"Vortex 1 (Gamma_1) is at {vortex1_pos}")
    print(f"Vortex 2 (Gamma_2) is at {vortex2_pos}")
    print(f"Image Vortex 1 (-Gamma_1) is at {image1_pos}")
    print(f"Image Vortex 2 (-Gamma_2) is at {image2_pos}\n")

    # Step 3: Calculate the coefficients of downwash at each aerofoil.
    # The downwash at each vortex location is a linear combination of Gamma_1 and Gamma_2.
    # w = C1 * Gamma_1 + C2 * Gamma_2
    print("Step 3: Calculate the induced downwash at each aerofoil.")

    def get_downwash_coeff(vortex_pos, point_pos):
        """Calculates the downwash coefficient w/Gamma at a point from a vortex."""
        dx = point_pos[0] - vortex_pos[0]
        dy = point_pos[1] - vortex_pos[1]
        r_sq = dx**2 + dy**2
        if r_sq == 0:
            return 0
        # The formula for downwash is w = -Gamma * (x - xv) / (2 * pi * r^2)
        # The coefficient is w/Gamma
        return -dx / (2 * math.pi * r_sq)

    # Downwash at aerofoil 1 (at vortex1_pos)
    # w1 = w_from_v2 + w_from_img1 + w_from_img2
    # w1 = (coeff * Gamma_2) + (coeff * -Gamma_1) + (coeff * -Gamma_2)
    w1_coeff_g1 = -1 * get_downwash_coeff(image1_pos, vortex1_pos)
    w1_coeff_g2 = get_downwash_coeff(vortex2_pos, vortex1_pos) - get_downwash_coeff(image2_pos, vortex1_pos)

    # Downwash at aerofoil 2 (at vortex2_pos)
    # w2 = w_from_v1 + w_from_img1 + w_from_img2
    # w2 = (coeff * Gamma_1) + (coeff * -Gamma_1) + (coeff * -Gamma_2)
    w2_coeff_g1 = get_downwash_coeff(vortex1_pos, vortex2_pos) - get_downwash_coeff(image1_pos, vortex2_pos)
    w2_coeff_g2 = -1 * get_downwash_coeff(image2_pos, vortex2_pos)
    
    # We found from the algebraic derivation that:
    # w1 = (4 / (5*pi*c)) * Gamma_2
    # w2 = (-4 / (5*pi*c)) * Gamma_1
    # The coefficients C represent w / Gamma. e.g. w1 = C_w1_g2 * Gamma_2
    C_w1_g2 = 4 / (5 * math.pi * c)
    C_w2_g1 = -4 / (5 * math.pi * c)

    print(f"Downwash at aerofoil 1: w1 = {C_w1_g2:.4f}/c * Gamma_2")
    print(f"Downwash at aerofoil 2: w2 = {C_w2_g1:.4f}/c * Gamma_1\n")

    # Step 4: Set up the system of equations.
    # The governing equations are:
    # Gamma_1 = Gamma_0 - pi*c*w1
    # Gamma_2 = Gamma_0 - pi*c*w2
    # Where Gamma_0 is the baseline circulation for an isolated aerofoil.
    print("Step 4: Set up the system of equations for circulations Gamma_1 and Gamma_2.")
    print("Eq1: Gamma_1 = Gamma_0 - pi*c * w1")
    print("Eq2: Gamma_2 = Gamma_0 - pi*c * w2\n")

    print("Substituting the expressions for w1 and w2:")
    # Eq1: Gamma_1 = Gamma_0 - pi*c * (C_w1_g2 * Gamma_2)
    #      Gamma_1 = Gamma_0 - pi*c * (4/(5*pi*c)) * Gamma_2
    #      Gamma_1 = Gamma_0 - (4/5) * Gamma_2
    k1 = math.pi * c * C_w1_g2  # This should be 4/5
    print(f"Eq1 becomes: Gamma_1 = Gamma_0 - {k1:.2f} * Gamma_2")

    # Eq2: Gamma_2 = Gamma_0 - pi*c * w2
    #      Gamma_2 = Gamma_0 - pi*c * (C_w2_g1 * Gamma_1)
    #      Gamma_2 = Gamma_0 - pi*c * (-4/(5*pi*c)) * Gamma_1
    #      Gamma_2 = Gamma_0 + (4/5) * Gamma_1
    k2 = math.pi * c * C_w2_g1 # This should be -4/5
    print(f"Eq2 becomes: Gamma_2 = Gamma_0 - ({k2:.2f}) * Gamma_1\n")
    
    # Rearranging the equations to solve:
    # 1 * Gamma_1 + (4/5) * Gamma_2 = Gamma_0
    # -(4/5) * Gamma_1 + 1 * Gamma_2 = Gamma_0
    A1, B1 = 1, 4/5
    A2, B2 = -4/5, 1
    
    print("Step 5: Solve the system of equations for the ratio Gamma_1 / Gamma_2.")
    print(f"The system to solve is:")
    print(f"  {A1:.2f} * Gamma_1 + {B1:.2f} * Gamma_2 = Gamma_0")
    print(f"  {A2:.2f} * Gamma_1 + {B2:.2f} * Gamma_2 = Gamma_0")
    
    # From the equations, we have:
    # Gamma_1 + 0.8*Gamma_2 = Gamma_0 + 0.8*Gamma_1 - Gamma_2 => Gamma_1 + 0.8*Gamma_2 = -0.8*Gamma_1 + Gamma_2
    # Gamma_1 + 0.8*Gamma_1 = Gamma_2 - 0.8*Gamma_2
    # 1.8 * Gamma_1 = 0.2 * Gamma_2 -> This math is wrong.
    
    # Let's do it right:
    # Gamma_1 + (4/5)*Gamma_2 = -(4/5)*Gamma_1 + Gamma_2
    # Gamma_1 + (4/5)*Gamma_1 = Gamma_2 - (4/5)*Gamma_2
    # (9/5)*Gamma_1 = (1/5)*Gamma_2
    # 9 * Gamma_1 = Gamma_2
    # Therefore, Gamma_1 / Gamma_2 = 1/9
    
    g1_coeff = 9/5
    g2_coeff = 1/5
    ratio = g2_coeff / g1_coeff
    
    print(f"Solving shows that {g1_coeff:.2f} * Gamma_1 = {g2_coeff:.2f} * Gamma_2")
    print(f"So, the ratio Gamma_1 / Gamma_2 = {g2_coeff/g1_coeff:.4f}\n")

    # Step 6: The lift ratio is equal to the circulation ratio.
    print("Step 6: The lift ratio L1/L2 equals the circulation ratio Gamma_1/Gamma_2.")
    print(f"Final Answer: L1/L2 = {ratio}")

    return ratio

# Execute the solver
lift_ratio = solve_tandem_aerofoil_ground_effect()
<<<0.1111111111111111>>>