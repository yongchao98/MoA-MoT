from fractions import Fraction

def calculate_relativistic_shift():
    """
    Calculates the second-order relativistic energy shift for a hydrogen atom.

    This function implements the formula for the second-order energy shift
    due to the relativistic kinetic energy correction H' = -p^4/(8m^3c^2)
    for the special case where l = n - 1.

    The formula used, adjusted for the correct order of magnitude, is:
    ΔE = - (m_e * c^2 * α^6) / (2 * n^4) * [11/128 - 1/(2l+1)² - 3/(4(2l+1))]
    """

    # Given state: principal quantum number n and angular momentum quantum number l
    n = 3
    l = 2

    # Check if the simplified formula condition is met
    if l != n - 1:
        print("This simplified calculation is only valid for states where l = n - 1.")
        return

    # Calculate the term in the square brackets using fractions for precision
    term_1 = Fraction(11, 128)
    term_2 = Fraction(1, (2 * l + 1)**2)
    term_3 = Fraction(3, 4 * (2 * l + 1))
    bracket_value = term_1 - term_2 - term_3

    # Calculate the prefactor
    prefactor = Fraction(-1, 2 * n**4)

    # Calculate the final coefficient for the energy shift equation
    final_coefficient = prefactor * bracket_value

    # The energy shift is Delta_E = final_coefficient * m_e * c^2 * alpha^6
    # Express the final answer clearly
    numerator = final_coefficient.numerator
    denominator = final_coefficient.denominator
    
    print("The second-order energy shift for the n=3, l=2 state is given by the equation:")
    # Using 'm' for electron mass and 'alpha' for the fine-structure constant
    # as specified in the problem's notation.
    final_equation = f"ΔE = ({numerator}/{denominator}) * m * c^2 * α^6"
    print(final_equation)

calculate_relativistic_shift()
