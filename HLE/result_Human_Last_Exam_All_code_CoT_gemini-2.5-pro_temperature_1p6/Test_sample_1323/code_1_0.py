import sympy as sp

def get_q1_expression():
    """
    Determines and prints the expression for the term '?_1'.
    
    This function uses the sympy library to create a symbolic representation
    of the mathematical result.
    """
    
    # Define symbolic variables
    # i and j are indices (1 or 2)
    # x1 and x2 are coordinates of the vector x
    i, j = sp.symbols('i j', integer=True, positive=True)
    x1, x2 = sp.symbols('x_1 x_2')
    
    # Define h as a general smooth function of x = (x1, x2)
    h = sp.Function('h')(x1, x2)
    
    # Define the Kronecker delta, which is 1 if i=j and 0 otherwise
    delta_ij = sp.KroneckerDelta(i, j)
    
    # Construct the expression for ?_1
    # ?_1 = (1/2) * delta_ij * h(x)
    one_half = sp.Rational(1, 2)
    q1 = one_half * delta_ij * h
    
    # Print the final result for ?_1
    print("The expression for the term ?_1 is:")
    
    # The final equation is ?_1 = ...
    # We output each component of the equation
    # The components are the symbol for ?_1, the equals sign, and the expression.
    q1_symbol = sp.Symbol('?_1')
    final_equation = sp.Eq(q1_symbol, q1)
    
    # Pretty print the equation
    sp.init_printing(use_unicode=True)
    print(sp.pretty(final_equation))
    
    # Also printing the raw expression to ensure all parts are visible
    # in any environment.
    print("\nRaw symbolic expression:")
    print(f"?_1 = {q1}")


if __name__ == '__main__':
    get_q1_expression()
