import sympy

def solve_tandem_aerofoil_ground_effect():
    """
    Calculates the lift ratio of two aerofoils in tandem formation
    in ground effect using the mirror image method.
    """
    # Define symbolic variables
    # Let c=1 for simplicity. U_inf and alpha are symbolic but will cancel out.
    c = 1.0
    U_inf, alpha = sympy.symbols('U_inf alpha')
    Gamma1, Gamma2 = sympy.symbols('Gamma1 Gamma2')

    # --- 1. Define Geometry ---
    s = 0.5 * c  # Separation
    h = 0.5 * c  # Ride height

    # Vortex and Control Point locations
    # Aerofoil 1 (front)
    x_v1 = 0.25 * c
    y_v1 = h
    x_cp1 = 0.75 * c
    y_cp1 = h

    # Aerofoil 2 (rear)
    x_v2 = c + s + 0.25 * c
    y_v2 = h
    x_cp2 = c + s + 0.75 * c
    y_cp2 = h
    
    # Image vortex locations
    x_iv1, y_iv1 = x_v1, -h
    x_iv2, y_iv2 = x_v2, -h
    
    print("--- Geometric Setup (assuming c=1) ---")
    print(f"Ride height h = {h}")
    print(f"Separation s = {s}")
    print(f"Vortex 1 at ({x_v1}, {y_v1}), CP1 at ({x_cp1}, {y_cp1})")
    print(f"Vortex 2 at ({x_v2}, {y_v2}), CP2 at ({x_cp2}, {y_cp2})")
    print(f"Image Vortex 1 at ({x_iv1}, {y_iv1})")
    print(f"Image Vortex 2 at ({x_iv2}, {y_iv2})\n")

    # --- 2. Induced Velocity Function ---
    # Calculates vertical velocity w induced by a vortex Gamma at (xv, yv) on a point (xp, yp)
    # The standard formula for downwash is w = -Gamma * (xp - xv) / (2 * pi * r^2)
    def induced_w(Gamma, xv, yv, xp, yp):
        dx = xp - xv
        dy = yp - yv
        r_sq = dx**2 + dy**2
        if r_sq == 0:
            return 0
        return -Gamma * dx / (2 * sympy.pi * r_sq)

    # --- 3. Setup Flow Tangency Equations ---
    # The sum of induced vertical velocities + freestream vertical velocity = 0
    # U_inf * alpha + sum(w_induced) = 0
    
    # Equation for Aerofoil 1 at CP1
    w_cp1_from_G2 = induced_w(Gamma2, x_v2, y_v2, x_cp1, y_cp1)
    w_cp1_from_iG1 = induced_w(-Gamma1, x_iv1, y_iv1, x_cp1, y_cp1) # from its own image
    w_cp1_from_iG2 = induced_w(-Gamma2, x_iv2, y_iv2, x_cp1, y_cp1) # from rear wing's image
    
    eq1 = sympy.Eq(U_inf * alpha + w_cp1_from_G2 + w_cp1_from_iG1 + w_cp1_from_iG2, 0)
    
    print("--- Equation for Aerofoil 1 ---")
    print("U_inf*alpha + w(Γ2->CP1) + w(-Γ1->CP1) + w(-Γ2->CP1) = 0")
    print(f"{U_inf*alpha} + {w_cp1_from_G2} + {w_cp1_from_iG1} + {w_cp1_from_iG2} = 0")
    print("Simplified Eq 1:", sympy.simplify(eq1), "\n")


    # Equation for Aerofoil 2 at CP2
    w_cp2_from_G1 = induced_w(Gamma1, x_v1, y_v1, x_cp2, y_cp2)
    w_cp2_from_iG1 = induced_w(-Gamma1, x_iv1, y_iv1, x_cp2, y_cp2) # from front wing's image
    w_cp2_from_iG2 = induced_w(-Gamma2, x_iv2, y_iv2, x_cp2, y_cp2) # from its own image

    eq2 = sympy.Eq(U_inf * alpha + w_cp2_from_G1 + w_cp2_from_iG1 + w_cp2_from_iG2, 0)

    print("--- Equation for Aerofoil 2 ---")
    print("U_inf*alpha + w(Γ1->CP2) + w(-Γ1->CP2) + w(-Γ2->CP2) = 0")
    print(f"{U_inf*alpha} + {w_cp2_from_G1} + {w_cp2_from_iG1} + {w_cp2_from_iG2} = 0")
    print("Simplified Eq 2:", sympy.simplify(eq2), "\n")

    # --- 4. Solve the System ---
    # We want the ratio Gamma1/Gamma2. We can solve for Gamma1 and Gamma2 in terms of U_inf*alpha
    # then take the ratio.
    solution = sympy.solve([eq1, eq2], [Gamma1, Gamma2])
    
    Gamma1_sol = solution[Gamma1]
    Gamma2_sol = solution[Gamma2]

    # The lift ratio L1/L2 is equal to the circulation ratio Gamma1/Gamma2
    lift_ratio = sympy.simplify(Gamma1_sol / Gamma2_sol)
    
    print("--- Solving the System ---")
    print("Solving for Γ1 and Γ2 and then finding the ratio L1/L2 = Γ1/Γ2.")
    print(f"Solution for Γ1 = {Gamma1_sol}")
    print(f"Solution for Γ2 = {Gamma2_sol}")
    print("\nThe final equation for the lift ratio is:")
    G1_num, G1_den = Gamma1_sol.as_numer_den()
    G2_num, G2_den = Gamma2_sol.as_numer_den()
    print(f"L1 / L2 = ({G1_num} / {G1_den}) / ({G2_num} / {G2_den}) = {lift_ratio}")
    
    return lift_ratio

final_ratio = solve_tandem_aerofoil_ground_effect()
print(f"\nFinal calculated lift ratio L1/L2 is: {float(final_ratio)}")

<<< -0.2 >>>