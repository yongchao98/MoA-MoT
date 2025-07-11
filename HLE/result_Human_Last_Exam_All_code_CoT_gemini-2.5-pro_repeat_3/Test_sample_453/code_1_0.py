import sympy

def calculate_relativistic_correction():
    """
    Calculates the first-order relativistic correction to the kinetic energy
    for a hydrogen atom in the state n=3, l=2.

    The first-order energy shift is given by:
    E_correction = <n,l,m| H' |n,l,m>
                 = - (1 / (2*m*c^2)) * < (E_n - V)^2 >
    where E_n is the unperturbed energy, and V is the Coulomb potential.
    This expands to:
    E_correction = - (1 / (2*m*c^2)) * [ (E_n)^2 - 2*E_n*<V> + <V^2> ]

    Using the Virial theorem for a 1/r potential, <V> = 2*E_n.
    The expectation value <V^2> can be related to E_n, n, and l.
    The final formula in terms of fundamental constants is:
    E_correction = - (m * c^2 * alpha^4 / (8 * n^4)) * (4*n / (l + 0.5) - 3)
    """

    # Given quantum numbers
    n = 3
    l = 2

    # We will calculate the numerical coefficient for the expression.
    # Let the expression be C * m * c^2 * alpha^4

    print(f"Calculating the first-order energy correction for n = {n}, l = {l}.")
    print("The formula for the energy shift is:")
    print("Delta_E = - (m * c^2 * alpha^4) / (8 * n^4) * ( (4*n / (l + 0.5)) - 3 )\n")

    # Calculate the term in the parenthesis
    term_in_parenthesis_val = (4 * n) / (l + 0.5) - 3
    print(f"The value of the term (4*n / (l + 0.5) - 3) is:")
    print(f"(4 * {n}) / ({l} + 0.5) - 3 = {term_in_parenthesis_val}\n")

    # Calculate the pre-factor coefficient
    pre_factor_denominator = 8 * (n**4)
    print(f"The value of the pre-factor -1 / (8 * n^4) is:")
    print(f"-1 / (8 * {n}^4) = -1 / {pre_factor_denominator}\n")

    # Calculate the total numerical coefficient
    # Using sympy to get a neat fraction
    total_coeff = sympy.Rational(-1, pre_factor_denominator) * sympy.Rational(str(term_in_parenthesis_val))
    
    # Extract numerator and denominator for the final print statement
    num = total_coeff.p
    den = total_coeff.q

    print("The total numerical coefficient is the product of these two parts:")
    print(f"Coefficient = (-1 / {pre_factor_denominator}) * {term_in_parenthesis_val} = {total_coeff}\n")

    # Final result expression
    print("The final expression for the energy shift is:")
    # The requirement is to output each number in the final equation.
    # The coefficient is -(1/360). So the numbers are 1 and 360.
    print(f"Delta_E = -({abs(num)}/{den}) * m * c^2 * alpha^4")

calculate_relativistic_correction()
<<<Delta_E = -(1/360) * m * c^2 * alpha^4>>>