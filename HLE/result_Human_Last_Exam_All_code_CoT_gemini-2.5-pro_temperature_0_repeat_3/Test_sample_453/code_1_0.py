import fractions

def calculate_relativistic_energy_shift():
    """
    Calculates the relativistic energy shift for a hydrogen atom in a given state.

    The problem asks for the energy shift for an electron in a hydrogen atom
    in a state with principal quantum number n=3 and angular momentum quantum
    number l=2. The perturbation is the relativistic correction to the kinetic
    energy.

    The energy shift is calculated using the first-order perturbation theory formula:
    Delta_E = - (m * c^2 * alpha^4) / (8 * n^4) * [ (4*n / (l + 1/2)) - 3 ]
    where m is the electron mass, c is the speed of light, and alpha is the
    fine-structure constant.
    """
    # Given quantum numbers
    n = 3
    l = 2

    # We will calculate the numerical coefficient of the expression.
    # The expression is of the form: Coefficient * m * c^2 * alpha^4

    # Calculate the term in the square brackets using fractions for precision
    term_in_brackets = (fractions.Fraction(4 * n) / (fractions.Fraction(l) + fractions.Fraction(1, 2))) - 3

    # Calculate the overall numerical coefficient
    # Coefficient = - (1 / (8 * n^4)) * term_in_brackets
    coefficient = -fractions.Fraction(1, 8 * n**4) * term_in_brackets

    # Extract numerator and denominator for the final output
    num = coefficient.numerator
    den = coefficient.denominator

    # Print the explanation and the final result
    print("Calculating the relativistic energy shift for a hydrogen atom.")
    print(f"State: n = {n}, l = {l}")
    print("The formula for the energy shift is:")
    print("ΔE = - (m*c²*α⁴)/(8n⁴) * [ 4n/(ℓ+1/2) - 3 ]")
    print("\nSubstituting the values:")
    print(f"The term in brackets is [ 4*{n}/({l}+0.5) - 3 ] = {float(term_in_brackets)}")
    print(f"The pre-factor is -1/(8*{n}⁴) = -1/{8*n**4}")
    print(f"The total numerical coefficient is (-1/{8*n**4}) * {term_in_brackets} = {coefficient}")

    print("\nThe final equation for the energy shift is:")
    # The final equation is printed showing each number
    print(f"ΔE = ({num}/{den}) * m * c² * α⁴")

calculate_relativistic_energy_shift()