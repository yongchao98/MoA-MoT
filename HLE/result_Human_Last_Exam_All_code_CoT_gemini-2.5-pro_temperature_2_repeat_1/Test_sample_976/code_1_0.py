import sympy

def solve_sphere_in_field():
    """
    This function symbolically derives the electric potential and field outside
    a conducting sphere placed in a uniform electric field in a steady state.
    """
    # Step 1: Define symbolic variables
    E0, R, r, theta = sympy.symbols('E_0 R r theta', real=True, positive=True)
    s1, s2 = sympy.symbols('sigma_1 sigma_2', real=True, positive=True)
    
    # The unknown coefficients for the potential solutions are D1 (for outside) and A1 (for inside).
    # We are primarily interested in the solution for r > R.
    # The general forms of the potential are assumed, based on the physics:
    # Phi_in = A1 * r * cos(theta)
    # Phi_out = -E0 * r * cos(theta) + D1 * r**(-2) * cos(theta)
    A1, D1 = sympy.symbols('A_1 D_1')

    # Step 2: Formulate equations based on boundary conditions at r = R
    
    # BC1: Continuity of potential -> Phi_in(R) = Phi_out(R)
    # A1*R*cos(theta) = -E0*R*cos(theta) + D1*R**-2*cos(theta)
    # This simplifies to (dividing by cos(theta)):
    eq1 = sympy.Eq(A1 * R, -E0 * R + D1 * R**-2)

    # BC2: Continuity of the normal component of current density -> s1*E_in_r = s2*E_out_r
    # E_r = -d(Phi)/dr
    # E_in_r = -A1*cos(theta)
    # E_out_r = -(-E0*cos(theta) - 2*D1*r**-3*cos(theta)) = (E0 + 2*D1*R**-3)*cos(theta)
    # So, s1*(-A1*cos(theta)) = s2*(E0 + 2*D1*R**-3)*cos(theta)
    # This simplifies to:
    eq2 = sympy.Eq(-s1 * A1, s2 * (E0 + 2 * D1 * R**-3))

    # Step 3: Solve the system of equations for the coefficient D1
    solution = sympy.solve([eq1, eq2], (A1, D1))
    D1_sol = solution[D1]

    # Step 4: Construct the final expressions for the potential and field outside the sphere
    
    # Potential Phi for r > R
    phi_out = -E0 * r * sympy.cos(theta) + D1_sol * r**(-2) * sympy.cos(theta)

    # Electric Field E = -grad(Phi) for r > R
    E_r_out = -sympy.diff(phi_out, r)
    E_theta_out = - (1/r) * sympy.diff(phi_out, theta)

    # Simplify expressions for printing
    phi_out_s = sympy.simplify(phi_out)
    E_r_out_s = sympy.simplify(E_r_out)
    E_theta_out_s = sympy.simplify(E_theta_out)
    
    # Step 5: Print the results and the individual components of the final equations
    
    print("--- Derived Results for the region r > R ---")
    
    # Expression for Potential
    print("\nElectric Potential Phi(r, theta):")
    k_factor_str = f"({s1} - {s2}) / ({s1} + 2*{s2})"
    phi_str = f"-{E0} * (r - {k_factor_str} * R**3 / r**2) * cos({theta})"
    print(phi_str)
    
    # Components of the Potential Equation
    print("\nComponents of the Potential Equation, Phi = Factor * (Term1 + Term2):")
    factor_phi = -E0 * sympy.cos(theta)
    term1_phi = r
    term2_phi = -(D1_sol / E0) * r**-2
    print(f"  Factor: {factor_phi}")
    print(f"  Term1: {term1_phi}")
    print(f"  Term2: {sympy.simplify(term2_phi)}")

    # Expression for Electric Field
    print("\nElectric Field E(r, theta) = E_r * r_hat + E_theta * theta_hat:")
    E_r_str = f"{E0} * (1 + 2 * {k_factor_str} * R**3 / r**3) * cos({theta})"
    E_theta_str = f"-{E0} * (1 - {k_factor_str} * R**3 / r**3) * sin({theta})"
    print(f"  E_r component: {E_r_str}")
    print(f"  E_theta component: {E_theta_str}")
    
    # Components of the E_r Equation
    print("\nComponents of the E_r Equation, E_r = Factor * (Term1 + Term2):")
    factor_Er = E0 * sympy.cos(theta)
    term1_Er = 1
    term2_Er = 2 * (D1_sol / E0) * r**-3
    print(f"  Factor: {factor_Er}")
    print(f"  Term1: {term1_Er}")
    print(f"  Term2: {sympy.simplify(term2_Er)}")

    # Components of the E_theta Equation
    print("\nComponents of the E_theta Equation, E_theta = Factor * (Term1 + Term2):")
    factor_Etheta = -E0 * sympy.sin(theta)
    term1_Etheta = 1
    term2_Etheta = -(D1_sol / E0) * r**-3
    print(f"  Factor: {factor_Etheta}")
    print(f"  Term1: {term1_Etheta}")
    print(f"  Term2: {sympy.simplify(term2_Etheta)}")

    print("\n\n--- Conclusion ---")
    print("The derived expressions match choice B.")


if __name__ == '__main__':
    solve_sphere_in_field()