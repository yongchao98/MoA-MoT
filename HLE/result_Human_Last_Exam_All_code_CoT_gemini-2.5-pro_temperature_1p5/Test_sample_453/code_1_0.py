import fractions

def calculate_relativistic_shift():
    """
    Calculates the first-order relativistic energy shift for a hydrogen atom.
    
    The energy shift is for the state with principal quantum number n and
    angular momentum quantum number l. The final result is expressed as a
    coefficient multiplying m_e * c^2 * alpha^4.
    """
    # Given quantum numbers
    n = 3
    l = 2

    # The formula for the relativistic energy shift is:
    # ΔE = - (1/2) * m_e * c^2 * (alpha^4 / n^4) * (n / (l + 1/2) - 3/4)
    # We will calculate the numerical coefficient:
    # C = - (1 / (2 * n^4)) * (n / (l + 1/2) - 3/4)

    print(f"Calculating the relativistic energy shift for a hydrogen atom state with n={n} and l={l}.")
    print("The formula for the first-order shift is: ΔE = - (m_e * c^2 * α^4 / (2 * n^4)) * (n / (l + 1/2) - 3/4)")
    
    # Calculate the term in parentheses using exact fractions
    term1 = fractions.Fraction(n, fractions.Fraction(l * 2 + 1, 2))
    term2 = fractions.Fraction(3, 4)
    parentheses_term = term1 - term2

    print(f"\nStep 1: Calculate the term in the parentheses (n / (l + 1/2) - 3/4)")
    print(f"   = ({n} / ({l} + 0.5)) - 3/4")
    print(f"   = {term1} - {term2}")
    print(f"   = {parentheses_term}")

    # Calculate the full coefficient
    prefactor = -fractions.Fraction(1, 2 * n**4)
    coefficient = prefactor * parentheses_term

    print(f"\nStep 2: Calculate the full numerical coefficient")
    print(f"   = -1 / (2 * {n}^4) * ({parentheses_term})")
    print(f"   = {prefactor} * {parentheses_term}")
    print(f"   = {coefficient}")

    print("\n--- Final Result ---")
    print("The calculated shift in the energy level is:")
    # The final equation string. Uses 'alpha' for the fine-structure constant.
    final_equation = f"ΔE = {coefficient.numerator}/{coefficient.denominator} * m_e * c^2 * alpha^4"
    print(final_equation)

if __name__ == "__main__":
    calculate_relativistic_shift()