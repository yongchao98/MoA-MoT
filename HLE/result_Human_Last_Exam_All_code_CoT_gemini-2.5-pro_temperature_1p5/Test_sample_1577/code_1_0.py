import sympy

def ground_space_degeneracy():
    """
    Calculates the symbolic ground space degeneracy for the toric code
    with n smooth holes and m rough holes.
    """
    # Define n and m as symbolic variables
    n = sympy.Symbol('n')
    m = sympy.Symbol('m')

    # The number of logical qubits, k, is (n-1) + (m-1)
    k = n + m - 2

    # The ground space degeneracy is 2^k
    # We will print this in a readable format.
    print("The ground space degeneracy (GSD) of the toric code is given by the formula:")
    
    # We want to print the equation GSD = 2^(n+m-2), including each number
    base = 2
    n_coeff = 1
    m_coeff = 1
    constant = -2
    
    # Print the equation part by part
    print(f"GSD = {base}^({n_coeff}*n + {m_coeff}*m + ({constant}))")
    
    # A cleaner print using sympy for the expression
    # This shows the logic more clearly.
    print("\nSymbolic representation of the exponent k:")
    print(f"k = (n - 1) + (m - 1) = {k}")

    print("\nFinal formula:")
    # Print the final expression 2^(n+m-2)
    final_expression_base = 2
    final_expression_exponent = k

    print(f"GSD = {final_expression_base}^({final_expression_exponent})")
    # To satisfy the output format of "each number in the final equation"
    print(f"\nFinal Equation with Numbers Explicit:")
    print(f"GSD = 2^(n + m - 2)")


if __name__ == "__main__":
    ground_space_degeneracy()
