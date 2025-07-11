import sympy

def calculate_energy_shift():
    """
    Calculates the first-order relativistic kinetic energy correction for a hydrogen atom.
    
    The formula for the energy shift is:
    ΔE = - (m_e * c^2 * α^4) / (2 * n^4) * (n / (l + 1/2) - 3/4)
    
    For the state n=3, l=2.
    """
    n = 3
    l = 2
    
    # Calculate the term in the parenthesis
    term = sympy.Rational(n, l + sympy.Rational(1, 2)) - sympy.Rational(3, 4)
    
    # Calculate the main coefficient
    coefficient = -sympy.Rational(1, 2 * n**4)
    
    # Total coefficient for m_e * c^2 * alpha^4
    total_coefficient = coefficient * term
    
    # Define symbols for the final expression
    m_e, c, alpha = sympy.symbols("m_e c alpha")
    
    # Create the symbolic expression for the energy shift
    energy_shift_expr = total_coefficient * m_e * c**2 * alpha**4
    
    # Get the numerator and denominator of the coefficient
    num, den = total_coefficient.p, total_coefficient.q
    
    # Print the final result step-by-step
    print(f"The energy shift ΔE for the state n={n}, l={l} is calculated as:")
    print(f"ΔE = - (m_e * c^2 * α^4) / (2 * n^4) * [n / (l + 1/2) - 3/4]")
    print(f"ΔE = - (m_e * c^2 * α^4) / (2 * {n}^4) * [{n} / ({l} + 1/2) - 3/4]")
    print(f"ΔE = - (m_e * c^2 * α^4) / (2 * {n**4}) * [{sympy.Rational(n, l + sympy.Rational(1,2))} - {sympy.Rational(3,4)}]")
    print(f"ΔE = - (m_e * c^2 * α^4) / ({2 * n**4}) * [{term}]")
    print(f"ΔE = ({coefficient}) * ({term}) * m_e * c^2 * α^4")
    print(f"The final expression for the energy shift is:")
    print(f"ΔE = ({num}/{den}) * m_e * c^2 * α^4")

calculate_energy_shift()
<<<ΔE = (-1/360) * m_e * c^2 * α^4>>>