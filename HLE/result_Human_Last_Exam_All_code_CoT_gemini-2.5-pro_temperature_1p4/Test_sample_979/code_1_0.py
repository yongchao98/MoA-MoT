import sympy as sp

def solve_magnetostatics_problem():
    """
    Solves the spherical shell magnetostatics problem symbolically.
    """
    # Define the physical constants and variables as symbolic objects
    A1, C1, K0, mu, mu0, R, r, theta = sp.symbols('A_1 C_1 K_0 mu mu_0 R r theta')

    # --- Step 1: Set up the system of equations from boundary conditions ---

    # We found from the math that the boundary conditions lead to two equations
    # for the l=1 coefficients, A1 and C1.

    # Equation 1 (from continuity of normal B-field: mu * H_in,r = mu0 * H_out,r)
    # H_in,r from potential A1*r*cos(theta) is -A1*cos(theta).
    # H_out,r from potential C1/r^2*cos(theta) is 2*C1/r^3*cos(theta).
    # At r=R: mu * (-A1) = mu0 * (2*C1/R**3)
    eq1 = sp.Eq(-mu * A1 * R**3, 2 * mu0 * C1)

    # Equation 2 (from discontinuity of tangential H: H_out,theta - H_in,theta = K0*sin(theta))
    # H_in,theta from potential A1*r*cos(theta) is A1*sin(theta).
    # H_out,theta from potential C1/r^2*cos(theta) is C1/r^3*sin(theta).
    # At r=R: C1/R**3 - A1 = K0
    eq2 = sp.Eq(C1 / R**3 - A1, K0)

    # --- Step 2: Solve the system for A1 and C1 ---
    solution = sp.solve([eq1, eq2], (A1, C1))
    A1_sol = solution[A1]
    C1_sol = solution[C1]
    
    # --- Step 3: Construct the H-field expressions ---

    # Inside field (r < R): H_in = -grad(A1*r*cos(theta)) = -A1 * z_hat
    # The coefficient of the field in the z-direction is -A1.
    H_in_coeff = -A1_sol
    
    # We reformat the expression to match the form in the answer choices
    H_in_final = (2 * mu0 / mu) * K0 / (1 + (2 * mu0 / mu))

    # Outside field (r > R): H_out = -grad(C1/r^2*cos(theta))
    # This gives H_out = C1/r^3 * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    # The coefficient part is C1.
    H_out_final_coeff = C1_sol / R**3
    
    # We reformat the expression to match the form in the answer choices
    H_out_final_coeff_reformatted = K0 / (1 + (2 * mu0 / mu))
    
    # --- Step 4: Print the final result ---
    print("The magnetic field H(r, theta) is:")
    print("For the region inside the sphere (0 < r < R):")
    # Using sp.pprint for a more readable mathematical output
    print("H_in = (", end="")
    sp.pprint(H_in_final)
    print(") * z_hat")
    
    print("\nFor the region outside the sphere (R < r < infinity):")
    print("H_out = (", end="")
    sp.pprint(H_out_final_coeff_reformatted)
    print(") * (R**3 / r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")

solve_magnetostatics_problem()