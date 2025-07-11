import sympy as sp

def solve_emi_shielding_force():
    """
    This function symbolically derives the force per unit area on the conducting plane
    using the magnetic scalar potential method and boundary conditions.
    """
    # Define symbols for the physical quantities and unknown coefficients
    K0, a, y, mu0, mu, x, d = sp.symbols('K_0 a y mu_0 mu x d', real=True, positive=True)
    C1, A2 = sp.symbols('C1 A2')

    # Define the magnetic scalar potential (psi_m) in the two regions
    # Region 1 (0 < x < d, air gap): The form satisfies H_x=0 at x=d.
    psi_m1 = C1 * sp.cosh(a * (x - d)) * sp.cos(a * y)
    
    # Region 2 (x < 0, magnetic material): The form satisfies H->0 as x->-infinity.
    psi_m2 = A2 * sp.exp(a * x) * sp.cos(a * y)

    # Calculate H = -grad(psi_m) for both regions
    H1x = -sp.diff(psi_m1, x)
    H1y = -sp.diff(psi_m1, y)
    H2x = -sp.diff(psi_m2, x)
    H2y = -sp.diff(psi_m2, y)

    # Set up equations from boundary conditions at x = 0
    # 1. Continuity of B_normal (mu0 * H1x = mu * H2x)
    # The cos(a*y) term is common and cancels out
    eq1_lhs = (mu0 * H1x / sp.cos(a * y)).subs(x, 0)
    eq1_rhs = (mu * H2x / sp.cos(a * y)).subs(x, 0)
    eq1 = sp.Eq(eq1_lhs, eq1_rhs)

    # 2. Discontinuity of H_tangential (H1y - H2y = K_z)
    # K_z = K0 * sin(a*y). The sin(a*y) term is common and cancels.
    eq2_lhs = ((H1y - H2y) / sp.sin(a * y)).subs(x, 0)
    eq2 = sp.Eq(eq2_lhs, K0)
    
    # Solve the system of two equations for the two unknown coefficients, C1 and A2
    solution = sp.solve([eq1, eq2], (C1, A2))
    C1_sol = solution[C1]
    
    # Calculate H_y at the conductor surface x = d
    H_y_at_conductor = H1y.subs(x, d)

    # Substitute the solved coefficient C1 into the expression for H_y
    H_y_final = H_y_at_conductor.subs(C1, C1_sol)
    
    # Calculate the force per unit area vector. Pressure is (1/2)*mu0*H^2.
    # Force is in the +x direction (pushing on the conductor).
    force_per_area_magnitude = (mu0 / 2) * H_y_final**2

    # Simplify the final expression
    force_per_area_simplified = sp.simplify(force_per_area_magnitude)

    # Extract numerator and denominator for clearer printing
    num, den = force_per_area_simplified.as_numer_denom()
    
    print("The derived force per unit area on the interface at x = d is:")
    print("Direction: +i_x")
    print("\n--- Components of the final equation ---")
    
    # We construct the printout to match the structure of the answer choices
    print(f"Coefficient: mu_0 / 2")
    
    # Isolate the parts of the formula for clarity
    K_term = K0**2 * sp.sin(a*y)**2
    denominator_term = (sp.cosh(a*d) + (mu0/mu)*sp.sinh(a*d))**2
    
    print(f"Numerator term: {K_term}")
    print(f"Denominator term: {denominator_term}")

    print("\n--- Final Expression ---")
    # Using sympy's pretty print for a nicely formatted output
    sp.pprint(force_per_area_simplified * sp.Symbol('i_x'))

solve_emi_shielding_force()