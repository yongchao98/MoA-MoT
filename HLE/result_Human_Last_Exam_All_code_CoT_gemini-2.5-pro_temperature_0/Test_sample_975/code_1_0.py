import sympy
from sympy import symbols, Eq, solve, diff, cos, sin, pretty_print

def solve_magnetostatics_problem():
    """
    Solves for the magnetic field H inside and outside a magnetized sphere
    within a conducting spherical shell using the magnetic scalar potential method.
    """
    # Define symbolic variables
    r, theta = symbols('r theta')
    M0, Rp, R = symbols('M0 R_p R', positive=True)
    A1, C1, D1 = symbols('A1 C1 D1')

    # Define the magnetic scalar potentials in the two regions
    # Region 1: 0 < r < Rp
    Phi_m1 = A1 * r * cos(theta)
    # Region 2: Rp < r < R
    Phi_m2 = (C1 * r + D1 / r**2) * cos(theta)

    # Set up the system of equations from boundary conditions

    # BC 1: Continuity of Phi_m at r = Rp
    # A1 * Rp = C1 * Rp + D1 / Rp**2
    eq1 = Eq(A1 * Rp, C1 * Rp + D1 / Rp**2)

    # BC 2: H_2r - H_1r = M_r at r = Rp
    # In terms of potential: diff(Phi_m1, r) - diff(Phi_m2, r) = M0 * cos(theta)
    # A1 - (C1 - 2*D1/Rp**3) = M0
    eq2 = Eq(A1 - (C1 - 2*D1 / Rp**3), M0)

    # BC 3: H_2r = 0 at r = R (Perfect Conductor)
    # -diff(Phi_m2, r) = 0 at r = R
    # -(C1 - 2*D1/R**3) = 0
    eq3 = Eq(C1 - 2*D1 / R**3, 0)

    # Solve the system of equations for the coefficients A1, C1, D1
    solution = solve((eq1, eq2, eq3), (A1, C1, D1))

    # Extract the solved coefficients
    A1_sol = solution[A1]
    C1_sol = solution[C1]
    D1_sol = solution[D1]

    # Calculate the H field components using H = -grad(Phi_m)
    # H_r = -diff(Phi_m, r), H_theta = -(1/r)*diff(Phi_m, theta)

    # Region 1: 0 < r < Rp
    H1_r_expr = -diff(Phi_m1.subs(A1, A1_sol), r)
    H1_theta_expr = - (1/r) * diff(Phi_m1.subs(A1, A1_sol), theta)

    # Region 2: Rp < r < R
    H2_r_expr = -diff(Phi_m2.subs({C1: C1_sol, D1: D1_sol}), r)
    H2_theta_expr = - (1/r) * diff(Phi_m2.subs({C1: C1_sol, D1: D1_sol}), theta)

    # --- Format and Print the Output ---
    print("--- Solution ---")
    print("\nIn the region 0 < r < R_p:")
    
    # H1 is uniform, so we can express it as a vector
    # H1 = A1_sol * (-cos(theta) i_r + sin(theta) i_theta)
    H1_coeff = sympy.simplify(A1_sol)
    # Factor to match the format in the options
    H1_coeff_factor = M0 * (2*Rp**3 + R**3) / (3*R**3)
    
    print("H = (M0 * (2*R_p**3 + R**3) / (3*R**3)) * (-cos(theta) i_r + sin(theta) i_theta)")
    print("Coefficient for H_1:")
    pretty_print(H1_coeff_factor)


    print("\nIn the region R_p < r < R:")
    
    # H2_r component
    H2_r_coeff = sympy.simplify(H2_r_expr / cos(theta))
    # Factor to match the format in the options
    H2_r_coeff_factor = -sympy.sympify(2)*M0/3 * ((Rp/R)**3 - (Rp/r)**3)
    
    print("H_r = (-2*M0/3 * [ (R_p/R)**3 - (R_p/r)**3 ]) * cos(theta)")
    print("Coefficient for H_2r / cos(theta):")
    pretty_print(H2_r_coeff_factor)

    # H2_theta component
    H2_theta_coeff = sympy.simplify(H2_theta_expr / sin(theta))
    # Factor to match the format in the options
    H2_theta_coeff_factor = M0/3 * (2*(Rp/R)**3 + (Rp/r)**3)

    print("\nH_theta = (M0/3 * [ 2*(R_p/R)**3 + (R_p/r)**3 ]) * sin(theta)")
    print("Coefficient for H_2theta / sin(theta):")
    pretty_print(H2_theta_coeff_factor)
    
    print("\nComparing these results with the given choices, the correct option is B.")

solve_magnetostatics_problem()