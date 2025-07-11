from fractions import Fraction

def calculate_relativistic_energy_shift():
    """
    Calculates the first-order relativistic kinetic energy correction for a
    hydrogen atom in the state (n, l).
    The formula for the energy shift is:
    E_shift = - (m * c^2 * alpha^4) / (2 * n^4) * (n / (l + 1/2) - 3/4)
    """
    # Given quantum numbers for the electron state
    n = 3
    l = 2
    
    print("This script calculates the first-order relativistic kinetic energy shift for an electron in a hydrogen atom.")
    print("The state is given by the principal quantum number n and the angular momentum quantum number l.")
    print(f"Given state: n = {n}, l = {l}\n")
    
    print("The formula for the energy shift E_shift is:")
    print("E_shift = - (m * c^2 * alpha^4) / (2 * n^4) * (n / (l + 0.5) - 3/4)\n")

    # Step-by-step calculation
    print("Step 1: Substitute the values of n and l into the formula.")
    print(f"E_shift = - (m * c^2 * alpha^4) / (2 * {n}^4) * ({n} / ({l} + 0.5) - 3/4)\n")

    # Evaluate the numerical parts
    n_fourth = n**4
    prefactor_den = 2 * n_fourth
    term1_num = n
    term1_den_float = l + 0.5
    term2 = Fraction(3, 4)
    
    term1 = Fraction(str(term1_num)) / Fraction(str(term1_den_float))
    
    print("Step 2: Evaluate the expression in the parenthesis.")
    print(f"Term 1: {n} / {l + 0.5} = {term1.numerator}/{term1.denominator}")
    print(f"Term 2: 3/4")
    
    parenthesis_result = term1 - term2
    print(f"Result of parenthesis: {term1.numerator}/{term1.denominator} - {term2.numerator}/{term2.denominator} = {parenthesis_result.numerator}/{parenthesis_result.denominator}\n")

    print("Step 3: Combine with the prefactor.")
    print(f"E_shift = - (m * c^2 * alpha^4) / {prefactor_den} * ({parenthesis_result.numerator}/{parenthesis_result.denominator})\n")

    # Final result
    total_coefficient = -Fraction(1, prefactor_den) * parenthesis_result
    num = total_coefficient.numerator
    den = total_coefficient.denominator

    print("Step 4: Calculate the final numerical coefficient.")
    print(f"E_shift = - (1 * {parenthesis_result.numerator}) / ({prefactor_den} * {parenthesis_result.denominator}) * (m * c^2 * alpha^4)")
    print(f"E_shift = - {parenthesis_result.numerator} / {prefactor_den * parenthesis_result.denominator} * (m * c^2 * alpha^4)")
    print(f"E_shift = {num}/{den} * (m * c^2 * alpha^4)\n")
    
    print("The final expression for the energy shift is:")
    print(f"{num}/{den} * m*c^2*alpha^4")
    
    return f"{num}/{den}"

calculate_relativistic_energy_shift()