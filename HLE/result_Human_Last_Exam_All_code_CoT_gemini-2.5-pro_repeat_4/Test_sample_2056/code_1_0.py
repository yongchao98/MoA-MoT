import sympy

def final_expression():
    """
    Constructs and prints the symbolic expression for l_k(n).
    """
    # Define the symbols n and k
    n = sympy.Symbol('n', real=True, positive=True)
    k = sympy.Symbol('k', real=True, positive=True)

    # Term 1: from the determinant of Sigma
    term1 = sympy.log(n + 1) / 2

    # Term 2: from the hyperbolic distance term
    term2 = (n - 1) * sympy.log(k)

    # Term 3: from the quadratic form
    term3 = -k**2 * (2 - 1/n)

    # The full expression for l_k(n)
    lk_n = term1 + term2 + term3

    # Print the final expression in a human-readable format
    # The numbers in the equation are: 1, 2, 1, 1, 2, 1
    print("The exact value of l_k(n) is:")
    sympy.pprint(lk_n, use_unicode=True)
    
    # As a single line string for clarity
    print("\nIn standard notation:")
    print("l_k(n) = (1/2)*ln(n + 1) + (n - 1)*ln(k) - k**2*(2 - 1/n)")


if __name__ == '__main__':
    final_expression()