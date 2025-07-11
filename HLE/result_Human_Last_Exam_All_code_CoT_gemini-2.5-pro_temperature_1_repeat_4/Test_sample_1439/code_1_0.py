import sympy as sp

def solve_critical_exponent_order():
    """
    This script symbolically determines the order in the coupling constant 'u'
    at which the critical exponent 'nu' first receives a non-vanishing contribution.
    """
    # 1. Define the symbolic variables.
    # 'u' is the coupling constant.
    # 'c' is a positive constant of proportionality from the one-loop calculation.
    u, c = sp.symbols('u c', positive=True)

    # 2. Define the relationship between nu and u.
    # In the epsilon-expansion, 1/nu = 2 - gamma_2(u).
    # To the lowest order in u, gamma_2(u) is linear in u: gamma_2(u) = c*u.
    # This is the one-loop contribution.
    gamma_2_u = c * u
    
    # The expression for nu as a function of u.
    nu_u = 1 / (2 - gamma_2_u)

    print("--- Analysis of the Critical Exponent ν(u) ---")
    print(f"The relationship between ν and the coupling constant u is given by: ν(u) = 1 / (2 - γ₂(u))")
    print(f"To the first-loop order, γ₂(u) = c*u. So, ν(u) = {nu_u}")
    print("-" * 50)

    # 3. Compute the Taylor series expansion of nu(u) around u=0.
    # We expand it to the second order in u to see the first and second contributions.
    series_nu = nu_u.series(u, 0, 3)

    print("The Taylor series expansion of ν(u) around u=0 is:")
    print(f"ν(u) ≈ {series_nu}")
    print("\nThe first term, ν(0), is the mean-field value.")
    
    mean_field_value = series_nu.subs(u, 0)
    print(f"Mean-field value ν(0) = {mean_field_value}")
    
    # 4. Identify the first non-vanishing contribution beyond the mean-field value.
    # The series terms are [1/2, c*u/4, c**2*u**2/8, O(u**3)]
    ordered_terms = series_nu.as_ordered_terms()
    
    if len(ordered_terms) > 2:
        # The first correction is the second term in the list.
        first_correction_term = ordered_terms[1]
        
        # Extract the power of 'u' from this term.
        # We can do this by creating a polynomial of the term and getting its degree.
        power_of_u = sp.Poly(first_correction_term, u).degree()

        print(f"\nThe first non-vanishing contribution (correction) to ν is the term: {first_correction_term}")
        
        # As requested, output the numbers in the equation for this term.
        # The term is of the form (c/A) * u^B
        coeff = first_correction_term / u
        denominator = sp.fraction(coeff)[1]

        print("\n--- Final Equation Details ---")
        print(f"This term can be written as (c/{denominator}) * u^{power_of_u}")
        print(f"The number in the denominator is: {denominator}")
        print(f"The power of the coupling constant 'u' is: {power_of_u}")
        print("-" * 50)

        print("\nTherefore, the critical exponent ν acquires its initial non-vanishing contribution at order:")
        print(power_of_u)
        
        # This is the final answer to be captured.
        return power_of_u
    else:
        print("Could not determine the order of the first correction.")
        return None

if __name__ == '__main__':
    solve_critical_exponent_order()
