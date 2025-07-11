import sympy

def find_max_n_expression():
    """
    This function determines the maximum value of n in terms of k
    based on the properties of the k-uniform intersecting family.

    The derivation shows that n_max = 2k - 1.
    This code uses the sympy library to represent this relationship symbolically.
    """

    # Define k as a symbolic variable to represent the size of the subsets
    k = sympy.symbols('k')

    # According to the mathematical proof, the maximum value for n is 2*k - 1.
    max_n_expr = 2 * k - 1

    # The final equation is n = 2*k - 1.
    # We print the components of this equation as requested.
    lhs = "n"
    op1 = "="
    num2 = 2
    op2 = "*"
    var_k = "k"
    op3 = "-"
    num1 = 1
    
    print("The maximum value of n in terms of k is given by the expression:")
    print(max_n_expr)
    print("\nThe final equation with each number and operator shown explicitly is:")
    print(f"{lhs} {op1} {num2} {op2} {var_k} {op3} {num1}")

# Execute the function to display the result.
find_max_n_expression()