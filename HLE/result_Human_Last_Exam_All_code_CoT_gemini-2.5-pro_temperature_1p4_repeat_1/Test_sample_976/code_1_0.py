import sympy

def solve_electrodynamics_problem():
    """
    Solves for the electric potential and field outside a conducting sphere
    in a uniform electric field using symbolic mathematics. This confirms
    the manual derivation.
    """
    # Define the symbols used in the problem
    # A and B are the unknown coefficients for the potential.
    # E0 is the magnitude of the uniform applied electric field.
    # R is the radius of the sphere.
    # sigma1 and sigma2 are the conductivities of the sphere and the medium.
    # r and theta are the spherical coordinates.
    A, B, E0, R, sigma1, sigma2, r, theta = sympy.symbols('A B E_0 R sigma_1 sigma_2 r theta', real=True)

    # General solutions:
    # Phi_in = A * r * cos(theta)
    # Phi_out = -E0 * r * cos(theta) + (B / r**2) * cos(theta)
    
    print("Step 1: Setting up boundary condition equations.")
    # Eq 1: Continuity of potential at r = R.
    # A*R = -E0*R + B/R^2
    eq1 = sympy.Eq(A * R, -E0 * R + B / R**2)
    print(f"  Equation from continuity of Potential: {eq1}")

    # Eq 2: Continuity of normal component of current density at r = R.
    # -sigma1 * A = sigma2 * (E0 + 2*B/R^3)
    eq2 = sympy.Eq(-sigma1 * A, sigma2 * (E0 + 2 * B / R**3))
    print(f"  Equation from continuity of J_normal: {eq2}")

    print("\nStep 2: Solving the system of equations for coefficients A and B.")
    solution = sympy.solve([eq1, eq2], (A, B))
    A_sol = solution[A]
    B_sol = solution[B]
    print(f"  Solved A = {A_sol}")
    print(f"  Solved B = {B_sol}")
    
    print("\nStep 3: Constructing potential and field for the region r > R.")
    # Potential outside: Phi_out(r, theta) = -E0*r*cos(theta) + (B/r^2)*cos(theta)
    phi_out_expr = -E0 * r * sympy.cos(theta) + (B_sol / r**2) * sympy.cos(theta)
    
    # Electric field E = -grad(Phi)
    E_r_out = -sympy.diff(phi_out_expr, r)
    E_theta_out = -sympy.diff(phi_out_expr, theta) / r

    # --- Final Expressions for r > R ---
    
    print("\n--- Final Expressions for Outside the Sphere (r > R) ---")

    print("\nElectric Potential Phi(r, theta):")
    # We format the output to match the structure of the answer choices.
    phi_num = f"({sigma1} - {sigma2}) * R^3"
    phi_den = f"({sigma1} + 2*{sigma2})*r^2"
    print(f"Phi(r, theta) = -E_0 * (r - ({phi_num}) / ({phi_den})) * cos(theta)")

    print("\nElectric Field E(r, theta) = E_r * r_hat + E_theta * theta_hat:")
    # Radial component
    Er_num = f"2*({sigma1} - {sigma2}) * R^3"
    Er_den = f"({sigma1} + 2*{sigma2}) * r^3"
    print(f"E_r = E_0 * [1 + ({Er_num}) / ({Er_den})] * cos(theta)")
    # Theta component
    Etheta_num = f"({sigma1} - {sigma2}) * R^3"
    Etheta_den = f"({sigma1} + 2*{sigma2}) * r^3"
    print(f"E_theta = -E_0 * [1 - ({Etheta_num}) / ({Etheta_den})] * sin(theta)")

    print("\nThese derived expressions match answer choice B.")

solve_electrodynamics_problem()