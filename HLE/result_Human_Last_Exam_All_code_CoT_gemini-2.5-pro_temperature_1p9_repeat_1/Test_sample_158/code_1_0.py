import sympy

def solve_magnetostatics_problem():
    """
    Solves for the required permeability of a cylindrical magnetic shell
    and the resulting interior magnetic field, ensuring no external field distortion.
    """
    # Define symbols used in the derivation.
    # H0: Magnitude of the external applied uniform magnetic field.
    # mu0: Permeability of free space.
    # R1, R2: Inner and outer radii of the cylindrical shell.
    # mu: Permeability of the magnetic material in the shell.
    # A, B, C: Coefficients for the scalar magnetic potential in the different regions.
    H0, mu0, R1, R2, mu = sympy.symbols('H_0 mu_0 R_1 R_2 mu', real=True, positive=True)
    A, B, C = sympy.symbols('A B C')
    
    # Based on the boundary conditions at the interfaces (rho=R1 and rho=R2),
    # we derive a system of four equations. Here, we directly use the
    # system of equations for the coefficients.

    # Equation 1 (from continuity of tangential H at rho=R2):
    eq1 = sympy.Eq(A + B/R2**2, -H0)
    
    # Equation 2 (from continuity of normal B at rho=R2):
    eq2 = sympy.Eq(mu * (A - B/R2**2), -mu0*H0)

    # Equation 3 (from continuity of tangential H at rho=R1):
    eq3 = sympy.Eq(C, A + B/R1**2)
    
    # Equation 4 (from continuity of normal B at rho=R1):
    eq4 = sympy.Eq(mu0*C, mu*(A - B/R1**2))
    
    # From equations 3 and 4, we eliminate C to get a constraint between A and B:
    # A + B/R1**2 = (mu/mu0)*(A - B/R1**2)
    # This simplifies to A*(mu - mu0) = B/R1**2 * (mu + mu0)
    constraint_eq = sympy.Eq(A * (mu - mu0), (B / R1**2) * (mu + mu0))
    
    # We solve equations 1 and 2 to find expressions for A and B.
    sol_AB = sympy.solve([eq1, eq2], [A, B])
    A_expr = sol_AB[A]
    B_expr = sol_AB[B]

    # Substitute these expressions for A and B into our constraint equation.
    final_eq_for_mu = constraint_eq.subs({A: A_expr, B: B_expr})
    
    # Simplify the equation to find the condition on mu.
    # The equation simplifies to H0*(mu**2 - mu0**2)*(R2**2 - R1**2) / (2*mu*R1**2) = 0.
    # For a non-trivial solution, the term (mu**2 - mu0**2) must be zero.
    # This gives mu = +mu0 or mu = -mu0.
    
    # The problem excludes the trivial case mu = mu0.
    # The non-trivial solution required is mu = -mu0.
    mu_solution_val = -mu0

    # Now, find the interior field H_int = -C * x_hat.
    # We need to calculate the coefficient C.
    C_expr = eq3.rhs # C = A + B/R1**2
    
    # Substitute the expressions for A and B into the equation for C.
    C_val = C_expr.subs({A: A_expr, B: B_expr})
    
    # Now substitute the required value for mu into the expression for C.
    C_val_final = C_val.subs(mu, mu_solution_val)
    C_val_final = sympy.simplify(C_val_final)
    
    # The interior field is H_int = -C * x_hat. Its magnitude is -C.
    H_int_magnitude = -C_val_final

    # --- Print the final results ---
    print("This script solves for the required material properties of a cylindrical shell to prevent distortion of an external magnetic field.")
    print("\n--- Problem Analysis ---")
    
    # We construct and display the key constraint equation derived from the physics.
    constraint_to_display = sympy.Eq((mu**2 - mu0**2) * (R2**2 / R1**2 - 1), 0)
    print("The physical boundary conditions lead to the following mathematical constraint:")
    print(f"Constraint Equation: {sympy.pretty(constraint_to_display, use_unicode=False)}")
    print("\nThis equation implies one of the following must be true:")
    print(f"1. mu = mu_0: The trivial case (excluded by the problem).")
    print(f"2. R1 = R2: A physically trivial case of a zero-thickness shell.")
    print(f"3. mu = -mu_0: The required non-trivial physical solution.")

    print("\n--- Final Answer ---")
    print("\n1. Required Permeability of the Shell Material")
    print("The required permeability `mu` for the shell material is:")
    print(f"mu = {sympy.pretty(mu_solution_val, use_unicode=False)}")
    print("In terms of relative permeability, mu_r = mu/mu_0:")
    print(f"mu_r = {sympy.pretty(mu_solution_val/mu0, use_unicode=False)}")

    print("\n2. Magnetic Field in the Interior Region")
    print("The magnetic field `H_int` in the interior region (rho < R1) is uniform, pointing in the x-direction:")
    # Create a symbolic representation of the vector for printing
    H_int_vector = f"({sympy.pretty(H_int_magnitude, use_unicode=False)}) * x_hat"
    print(f"H_int = {H_int_vector}")

if __name__ == '__main__':
    solve_magnetostatics_problem()
    # The expression for mu is `mu = -mu_0`.
    # The expression for H_int is `H0*(R2/R1)**2 * x_hat`.
    # Let's write the final values in the required format for the submission.
    # The solution for mu in terms of mu0: mu_r = -1.
    # The solution for H_int in terms of H0, R1, R2: H_int = H0*(R2/R1)**2 * x_hat
    
    # Based on this, we extract the core answers requested by the prompt.
    mu_r_answer = -1
    # H_int_answer depends on H0, R1, R2, so it's a symbolic expression.
    # The prompt will be answered by the script's printout.
    # No specific numerical values are given for R1, R2, H0 so this is the final form.
    # If forced to pick a single value, it could be the relative permeability.
    # >>> mu_r_answer -> -1
    # <<< mu_r_answer = -1 ; H_int = H0*(R2/R1)**2 * x_hat >>>
    # The output format requests one answer. The question asks for two things.
    # Let's provide both in a structured way.
    mu = sympy.Symbol('mu')
    mu0 = sympy.Symbol('mu_0')
    H0 = sympy.Symbol('H_0')
    R1 = sympy.Symbol('R_1')
    R2 = sympy.Symbol('R_2')
    H_int = H0 * (R2/R1)**2
    permeability_relation = sympy.Eq(mu, -mu0)
    final_answer_string = f"Required permeability: {permeability_relation}; Interior field magnitude: H_int_mag/H_0 = (R2/R1)**2"
    # To satisfy the simple answer format, providing mu_r is probably the best single numeric answer.
    # <<<-1>>> is probably too simplistic. Let's provide the final symbolic results as structured text.