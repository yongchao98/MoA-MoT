import sympy

def find_term_q1():
    """
    This script represents and prints the mathematical expression for the term '?_1'.
    The term is derived from the second partial derivatives of the 2D potential integral.
    The derivation shows that ?_1 = (1/2) * h(x) * delta_ij.
    This script uses the sympy library to represent this symbolic result.
    """

    # Define the necessary mathematical symbols
    h = sympy.Function('h')
    x1, x2 = sympy.symbols('x1 x2')
    i, j = sympy.symbols('i j', integer=True)

    # The derived coefficient is 1/2.
    # We define it as a sympy rational number to preserve its form.
    coefficient = sympy.S(1)/2
    
    # The term involves the function h(x)
    h_of_x = h(x1, x2)
    
    # The term involves the Kronecker delta, which is 1 if i=j and 0 otherwise.
    kronecker_delta = sympy.KroneckerDelta(i, j)

    # The complete expression for ?_1
    question_mark_1 = coefficient * h_of_x * kronecker_delta

    # Print the final expression in a formatted way.
    # The numbers in the fraction 1/2 are explicitly printed.
    print("The term ?_1 is determined to be:")
    print(f"({sympy.numer(coefficient)}/{sympy.denom(coefficient)}) * h({x1}, {x2}) * delta({i}, {j})")

    print("\nWhere:")
    print(f" - h({x1}, {x2}) is the function h evaluated at the point x = ({x1}, {x2}).")
    print(f" - delta({i}, {j}) is the Kronecker delta.")

# Run the function to display the result
find_term_q1()