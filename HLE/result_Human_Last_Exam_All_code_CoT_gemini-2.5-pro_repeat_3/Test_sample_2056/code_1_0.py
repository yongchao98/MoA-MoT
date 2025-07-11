import sympy as sp

def solve_and_print_formula():
    """
    This function derives and prints the symbolic formula for l_k(n).
    The derivation steps are outlined in the text above.
    """
    n = sp.Symbol('n', integer=True, positive=True)
    k = sp.Symbol('k', real=True, positive=True)

    # The derived formula for l_k(n)
    # l_k(n) = (1/2)*ln(n+1) - k^2*(2 - 1/n) + (n-1)*ln(k)
    term1_coeff_num = 1
    term1_coeff_den = 2
    term2_const = 2
    term2_coeff_num = 1
    term3_const = 1

    term1 = sp.Rational(term1_coeff_num, term1_coeff_den) * sp.log(n + 1)
    term2 = -k**2 * (term2_const - sp.Rational(term2_coeff_num, n))
    term3 = (n - term3_const) * sp.log(k)

    final_expression = term1 + term2 + term3

    print("The exact value of l_k(n) in terms of n and k is:")
    sp.pprint(final_expression)

    # As requested by the prompt, outputting the numbers in the final equation
    print("\nThe numbers in the final equation are:")
    print(f"Coefficient of ln(n+1): {term1_coeff_num}/{term1_coeff_den}")
    print(f"Constant in k^2 term: {term2_const}")
    print(f"Coefficient of 1/n in k^2 term: {term2_coeff_num}")
    print(f"Constant subtracted from n in ln(k) term: {term3_const}")


if __name__ == '__main__':
    solve_and_print_formula()
