import sympy

def solve_for_q1():
    """
    This function determines and prints the expression for the term ?_1
    based on the mathematical derivation.
    """
    # Define symbolic variables and the function h(x)
    x_1, x_2 = sympy.symbols('x_1 x_2')
    x = sympy.Matrix([x_1, x_2])
    h = sympy.Function('h')(x)

    # Define symbolic indices i and j
    i, j = sympy.symbols('i j', integer=True, positive=True, constants=True)

    # The constant factor in the expression for ?_1
    C = sympy.Rational(1, 2)
    
    # The Kronecker delta
    delta_ij = sympy.KroneckerDelta(i, j)
    
    # Construct the expression for ?_1
    q1 = C * h * delta_ij
    
    # Print the general result for ?_1
    print("The term ?_1 is given by the following expression:")
    sympy.pprint(q1, use_unicode=True)
    
    print("\nThis means:")
    # Display the result for the case i = j
    print("If i = j:")
    q1_equal = q1.subs(j, i)
    sympy.pprint(q1_equal, use_unicode=True)

    # Display the result for the case i != j
    print("\nIf i != j:")
    # Sympy automatically evaluates KroneckerDelta(i,j) to 0 if it can prove i != j.
    # We create a new delta where indices are unequal to show the result.
    q1_unequal = C * h * sympy.KroneckerDelta(1, 2)
    sympy.pprint(q1_unequal, use_unicode=True)

if __name__ == '__main__':
    solve_for_q1()