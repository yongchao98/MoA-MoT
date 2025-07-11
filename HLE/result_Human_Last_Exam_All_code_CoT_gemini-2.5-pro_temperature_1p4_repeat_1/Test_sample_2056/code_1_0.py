import sympy

def solve_lk_n():
    """
    This function calculates the symbolic expression for l_k(n) and prints its components.
    """
    # Define the symbols used in the formula
    n = sympy.Symbol('n', integer=True, positive=True)
    k = sympy.Symbol('k', real=True, positive=True)

    # The final derived formula for l_k(n)
    # l_k(n) = (1/2)*ln(n+1) + (n-1)*ln(k) - (1/2)*ln(2) - k**2 * (2 - 1/n)
    
    # Let's build the expression term by term to identify the numbers
    term1_coeff = sympy.Rational(1, 2)
    term1_arg = n + 1
    
    term2_coeff_n = n - 1

    term3_coeff = -sympy.Rational(1, 2)
    term3_arg = 2
    
    term4_factor_k = -k**2
    term4_arg_n = 2 - sympy.Rational(1, n)

    # Combine the terms to form the final expression
    l_k_n_expression = term1_coeff * sympy.log(term1_arg) + \
                       term2_coeff_n * sympy.log(k) + \
                       term3_coeff * sympy.log(term3_arg) + \
                       term4_factor_k * term4_arg_n

    print("The exact value of l_k(n) in terms of n and k is:")
    sympy.pprint(l_k_n_expression)

    # Output each number in the final equation as requested by the prompt
    # The structure is A*ln(n+B) + (n-C)*ln(k) - D*ln(E) - k^2*(F - G/n)
    A = term1_coeff
    B = 1
    C = 1
    D = -term3_coeff
    E = term3_arg
    F = 2
    G = 1

    print("\n--- Equation Breakdown ---")
    print(f"The equation has the form: A*ln(n + B) + (n - C)*ln(k) - D*ln(E) - k^2*(F - G/n)")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print(f"F = {F}")
    print(f"G = {G}")


if __name__ == '__main__':
    solve_lk_n()
