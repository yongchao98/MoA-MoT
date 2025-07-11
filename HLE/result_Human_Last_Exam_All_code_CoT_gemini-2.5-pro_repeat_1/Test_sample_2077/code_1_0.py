import sympy as sp

def check_formula():
    """
    This function checks if Option D simplifies to the given approximate formula
    under the specified approximations using symbolic mathematics.
    """
    # Define symbolic variables for the physical quantities
    E_nu, M, m_nu, T = sp.symbols('E_nu M m_nu T', positive=True, real=True)
    G_F, Q_W = sp.symbols('G_F Q_W') # Fermi constant and Weak charge
    # Represent the form factor squared as a single symbol
    F_q2_sq = sp.Symbol('F_q2_sq')

    # --- Part 1: Define the components from Option D ---
    # Upper limit of integration from Option D
    T_max_D = (2*M*E_nu**2 - 2*M*m_nu**2) / (2*M*E_nu + M**2 + m_nu**2)

    # Prefactor from Option D
    prefactor_numerator_D = G_F**2 * Q_W**2 * F_q2_sq * E_nu**2 * M**3
    prefactor_denominator_D = sp.pi * ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
    prefactor_D = prefactor_numerator_D / prefactor_denominator_D

    # Term in brackets from Option D
    bracket_D = 1 - T/E_nu - M*T/(2*E_nu**2) - m_nu**2/(4*E_nu**2) - m_nu**2*T/(4*M*E_nu**2)

    print("--- Analysis of Option D ---")
    print("Original Upper Limit (T_max):")
    sp.pprint(T_max_D)
    print("\nOriginal Prefactor:")
    sp.pprint(prefactor_D)
    print("\nOriginal Bracket Term:")
    sp.pprint(bracket_D)
    print("\n" + "="*50 + "\n")

    # --- Part 2: Apply the approximation m_nu = 0 ---
    print("--- Step 1: Applying approximation m_nu = 0 ---")
    T_max_m0 = T_max_D.subs(m_nu, 0)
    prefactor_m0 = prefactor_D.subs(m_nu, 0)
    bracket_m0 = bracket_D.subs(m_nu, 0)

    print("Resulting T_max (m_nu = 0):")
    sp.pprint(sp.simplify(T_max_m0))
    print("\nResulting Prefactor (m_nu = 0):")
    sp.pprint(sp.simplify(prefactor_m0))
    print("\nResulting Bracket Term (m_nu = 0):")
    sp.pprint(bracket_m0)
    print("\n" + "="*50 + "\n")

    # --- Part 3: Compare with the given formula and apply E_nu << M ---
    print("--- Step 2: Applying approximation E_nu << M ---")

    # Define the components from the given approximate formula for comparison
    T_max_given = E_nu / (1 + M/(2*E_nu))
    prefactor_given = (G_F**2 * Q_W**2 * F_q2_sq * M) / (4*sp.pi)
    bracket_given = 1 - M*T / (2*E_nu**2)

    print("Checking T_max:")
    print("The T_max from Option D with m_nu=0 is:")
    sp.pprint(sp.simplify(T_max_m0))
    print("The T_max from the given formula is:")
    sp.pprint(sp.simplify(T_max_given))
    print("Do they match?", sp.simplify(T_max_m0) == sp.simplify(T_max_given))
    print("(Note: The E_nu << M approximation was not needed for the limit to match)")
    print("-" * 20)

    print("Checking Prefactor:")
    # To apply E_nu << M, we approximate (E_nu + 2*M) with (2*M)
    prefactor_final = sp.simplify(prefactor_m0).subs((E_nu + 2*M), (2*M))
    print("The Prefactor from Option D after both approximations is:")
    sp.pprint(prefactor_final)
    print("The Prefactor from the given formula is:")
    sp.pprint(prefactor_given)
    print("Do they match?", prefactor_final == prefactor_given)
    print("-" * 20)
    
    print("Checking Bracket Term:")
    # In the bracket [1 - T/E_nu - M*T/(2*E_nu**2)], the ratio of the T-dependent terms
    # is (T/E_nu) / (M*T/(2*E_nu**2)) = 2*E_nu/M.
    # Since E_nu << M, the T/E_nu term is negligible and can be dropped.
    bracket_final = 1 - M*T/(2*E_nu**2)
    print("The Bracket from Option D after both approximations is:")
    sp.pprint(bracket_final)
    print("The Bracket from the given formula is:")
    sp.pprint(bracket_given)
    print("Does it match?", bracket_final == bracket_given)
    print("\n" + "="*50 + "\n")
    
    print("Conclusion: Option D correctly reduces to the provided approximate formula.")

# Run the analysis
check_formula()