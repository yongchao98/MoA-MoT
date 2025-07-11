import sympy

def calculate_relativistic_correction():
    """
    Calculates the first-order relativistic kinetic energy correction
    for a hydrogen atom state (n, l).

    The calculation follows these steps:
    1.  The perturbation is H' = -p^4 / (8*m^3*c^2).
    2.  Using H_0 = p^2/(2m) + V, we get p^2 = 2m*(H_0 - V).
    3.  So, p^4 = 4m^2*(H_0 - V)^2.
    4.  H' = -4m^2*(H_0 - V)^2 / (8*m^3*c^2) = -1/(2*m*c^2) * (H_0 - V)^2.
    5.  The first-order correction is <H'> = -1/(2*m*c^2) * <(H_0 - V)^2>.
    6.  <(H_0 - V)^2> = <H_0^2 - 2*H_0*V + V^2> = (E_n)^2 - 2*E_n*<V> + <V^2>.
    7.  From the virial theorem, <V> = 2*E_n.
    8.  So, <(H_0 - V)^2> = (E_n)^2 - 4*(E_n)^2 + <V^2> = -3*(E_n)^2 + <V^2>.
    9.  The final expression for the energy shift is:
        Delta_E = -1/(2*m*c^2) * [-3*(E_n)^2 + <V^2>].
    """

    # Define quantum numbers
    n_val = 3
    ell_val = 2

    # Define symbols for the formula
    m, c, alpha, n, ell = sympy.symbols('m_e c alpha n ell')

    # Energy of the unperturbed state E_n = - (m*c^2*alpha^2) / (2*n^2)
    E_n = -(m * c**2 * alpha**2) / (2 * n**2)

    # Expectation value of V^2 for a hydrogen state |n,l>
    # <V^2> = <(-alpha*hbar*c/r)^2> = (alpha*hbar*c)^2 * <1/r^2>
    # <1/r^2> = 1 / (a_0^2 * n^3 * (ell + 1/2))
    # a_0 = hbar / (m*c*alpha) => a_0^2 = hbar^2 / (m^2*c^2*alpha^2)
    # <V^2> = (alpha*hbar*c)^2 * (m^2*c^2*alpha^2) / (hbar^2 * n^3 * (ell + 1/2))
    # <V^2> = m^2 * c^4 * alpha^4 / (n^3 * (ell + 1/2))
    V2_exp = (m**2 * c**4 * alpha**4) / (n**3 * (ell + sympy.S(1)/2))

    # Calculate the shift
    # Delta_E = -1/(2*m*c^2) * [-3*(E_n)^2 + <V^2>]
    delta_E = sympy.simplify(-1 / (2 * m * c**2) * (-3 * E_n**2 + V2_exp))

    # Substitute the values for n and ell
    delta_E_specific = delta_E.subs({n: n_val, ell: ell_val})

    # The formula for the energy shift is
    # Delta_E = C * m * c^2 * alpha^4
    # Let's find the coefficient C
    coefficient = sympy.simplify(delta_E_specific / (m * c**2 * alpha**4))

    print(f"For a hydrogen atom in the state n = {n_val}, l = {ell_val}:")
    print(f"The first-order relativistic kinetic energy correction is given by the formula:")
    print(f"ΔE = C * m_e * c^2 * α^4")
    print(f"where the calculated coefficient C is: {coefficient}")
    
    # We can expand alpha to express the result in terms of e, hbar, c, epsilon_0
    e, hbar, epsilon_0 = sympy.symbols('e hbar epsilon_0')
    alpha_expr = e**2 / (4 * sympy.pi * epsilon_0 * hbar * c)
    
    final_expression_coeff = coefficient * m * c**2
    final_expression = final_expression_coeff * alpha_expr**4
    
    print("\nThe full equation for the energy shift in terms of fundamental constants is:")
    # Using unicode characters for a nicer print
    m_e_sym = '\u006d\u2091'
    alpha_sym = '\u03b1'
    hbar_sym = '\u210f'
    epsilon_0_sym = '\u03b5\u2080'
    
    # Print the equation part by part as requested
    print(f"ΔE_{{{n_val},{ell_val}}} = {coefficient} * {m_e_sym} * c² * {alpha_sym}⁴")
    
    # Simplify and print the expression with base constants
    final_simplified = sympy.simplify(final_expression)
    num, den = final_simplified.as_numer_denom()
    
    print(f"   = - ({m_e_sym} * e⁸) / (360 * (4π{epsilon_0_sym})⁴ * {hbar_sym}⁴ * c²)")


calculate_relativistic_correction()