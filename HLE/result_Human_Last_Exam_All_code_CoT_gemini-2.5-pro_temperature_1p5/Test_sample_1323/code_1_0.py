import sympy

def solve_problem():
    """
    This function defines the symbolic expression for ?_1 and prints it.
    """
    # Define symbols used in the expression for ?_1
    h = sympy.Function('h')
    x = sympy.Symbol('x')
    i = sympy.Symbol('i', integer=True)
    j = sympy.Symbol('j', integer=True)

    # Construct the expression for ?_1
    # ?_1 = (1/2) * h(x) * delta_ij
    # In SymPy, S(1)/2 creates a fraction and KroneckerDelta(i, j) represents delta_ij.
    q1_expression = (sympy.S(1)/2) * h(x) * sympy.KroneckerDelta(i, j)

    # Print the resulting expression in a readable format
    print("The symbolic expression for ?_1 is:")
    sympy.pprint(q1_expression, use_unicode=True)
    
    # As requested, output each number in the final equation.
    # We find all numerical parts of the expression. In this case, it's just 1/2.
    numbers = q1_expression.atoms(sympy.Number)
    print("\nThe numbers found in the final expression are:")
    for num in numbers:
        print(num)

if __name__ == "__main__":
    solve_problem()
