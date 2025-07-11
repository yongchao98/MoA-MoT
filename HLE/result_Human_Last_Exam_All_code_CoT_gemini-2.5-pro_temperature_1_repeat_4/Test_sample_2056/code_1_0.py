import sympy

def calculate_exact_value():
    """
    This function prints the final derived symbolic expression for l_k(n).
    The derivation is explained in the text above.
    """
    
    # Define symbolic variables for the expression
    n, k = sympy.symbols('n k')

    # The final derived formula for l_k(n) is:
    # l_k(n) = (1/2)ln(n+1) + (n-1)ln(k) - k^2(2 - 1/n)
    
    # To satisfy the prompt "output each number in the final equation",
    # we will construct the final output string manually.
    
    equation_part_1 = "(1/2) * log(n + 1)"
    equation_part_2 = "(n - 1) * log(k)"
    equation_part_3 = "k**2 * (2 - 1/n)"

    final_equation_string = f"l_k(n) = {equation_part_1} + {equation_part_2} - {equation_part_3}"

    print("The exact value of l_k(n) in terms of n and k is:")
    print(final_equation_string)

    # For verification and a more standard mathematical output, we can use sympy's pretty print.
    # print("\nSymbolic expression using sympy:")
    # final_expression = sympy.Rational(1, 2) * sympy.log(n + 1) + (n - 1) * sympy.log(k) - k**2 * (2 - 1/n)
    # sympy.pprint(final_expression)


if __name__ == '__main__':
    calculate_exact_value()