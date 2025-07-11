import sympy

def solve_electostatics_problem():
    """
    This function solves for the electric potential and field outside a conducting
    sphere in a uniform electric field using symbolic mathematics.
    """
    # Step 1: Define symbolic variables
    r, R, theta = sympy.symbols('r R theta', real=True, positive=True)
    E0, s1, s2 = sympy.symbols('E_0 sigma_1 sigma_2', real=True, positive=True)
    A1, D1 = sympy.symbols('A1 D1') # Unknown coefficients

    # Cos(theta) is the first Legendre Polynomial P1
    P1 = sympy.cos(theta)

    # Step 2: Define the general forms of the potential
    # Inside the sphere (r < R), potential must be finite at r=0
    # Phi_in = A1 * r * P1 (we only need the l=1 term due to the external field's form)
    Phi_in = A1 * r * P1

    # Outside the sphere (r > R), potential approaches -E0*r*cos(theta) at infinity
    # Phi_out = -E0*r*P1 + D1/r^2 * P1
    Phi_out = -E0 * r * P1 + D1 / (r**2) * P1

    # Step 3: Apply boundary conditions at r = R
    # Condition 1: Potential is continuous at r = R
    # Phi_in(R) = Phi_out(R)
    eq1 = sympy.Eq(Phi_in.subs(r, R), Phi_out.subs(r, R))

    # Condition 2: Normal component of current density is continuous at r = R
    # J_n = sigma * E_r = -sigma * d(Phi)/dr
    # J_1n(R) = J_2n(R)  => s1 * d(Phi_in)/dr |_(r=R) = s2 * d(Phi_out)/dr |_(r=R)
    dPhi_in_dr = sympy.diff(Phi_in, r)
    dPhi_out_dr = sympy.diff(Phi_out, r)
    eq2 = sympy.Eq(s1 * dPhi_in_dr.subs(r, R), s2 * dPhi_out_dr.subs(r, R))

    # Step 4: Solve the system of equations for the coefficients A1 and D1
    solution = sympy.solve([eq1, eq2], (A1, D1))
    
    # Extract the solved coefficients
    D1_sol = solution[D1]

    # Step 5: Substitute the solved coefficient back into the expression for Phi_out
    Phi_out_final = Phi_out.subs(D1, D1_sol)

    # Step 6: Calculate the electric field outside the sphere, E = -grad(Phi)
    # E_r = -d(Phi)/dr
    # E_theta = -(1/r) * d(Phi)/d(theta)
    Er_out = -sympy.diff(Phi_out_final, r)
    Etheta_out = -sympy.simplify(1/r * sympy.diff(Phi_out_final, theta))

    # Step 7: Print the results for the region outside the sphere (r > R)
    print("This problem is analogous to a dielectric sphere in a uniform electric field, with the substitution epsilon -> sigma.")
    print("The governing equation is Laplace's equation, solved with boundary conditions at r=R.")
    print("\n--- Results for the region outside the sphere (r > R) ---")
    
    # Format the potential for readability
    dipole_term = D1_sol / E0
    Phi_out_display = -E0 * (r - dipole_term.simplify() / r**2) * sympy.cos(theta)
    
    print("\nElectric Potential Phi(r, theta):")
    sympy.pprint(Phi_out_display, use_unicode=False)
    # This matches the form: -E_0 * (r - (sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**2)) * cos(theta)

    print("\nElectric Field E(r, theta) = E_r * r_hat + E_theta * theta_hat:")
    print("\nRadial component (E_r):")
    sympy.pprint(Er_out.simplify(), use_unicode=False)
    # This matches the form: E_0 * (1 + 2*(sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**3)) * cos(theta)
    
    print("\nAngular component (E_theta):")
    sympy.pprint(Etheta_out, use_unicode=False)
    # This matches the form: -E_0 * (1 - (sigma_1 - sigma_2)*R**3 / ((sigma_1 + 2*sigma_2)*r**3)) * sin(theta)

    print("\nComparing these results with the given options, we find that option B is the correct choice.")

solve_electostatics_problem()