import sympy

def solve_factorization():
    """
    This function calculates and prints a left coprime factorization for the given transfer function H(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # The transfer function matrix H(s)
    # H = sympy.Matrix([[(s-1)/(s+1), 1], [2/(s**2-1), 0]])

    # A left coprime factorization is H(s) = D(s)^-1 * N(s).
    # We found the following polynomial matrices D(s) and N(s).
    
    # Define the matrix D(s)
    D = sympy.Matrix([[s + 1, 0],
                      [1, s - 1]])

    # Define the matrix N(s)
    N = sympy.Matrix([[s - 1, s + 1],
                      [1, 1]])

    # Print the results in a structured format
    print("A left coprime factorization of the transfer function H(s) is given by:")
    print("H(s) = D(s)^-1 * N(s)")
    
    print("\nwhere D(s) is:")
    sympy.pprint(D, use_unicode=False)
    
    print("\nand N(s) is:")
    sympy.pprint(N, use_unicode=False)
    
    # As requested, outputting each polynomial element of the matrices
    print("\n--- Details of each element in the final equation ---")
    
    print("\nMatrix D(s):")
    print(f"D_11(s) = {D[0,0]}")
    print(f"D_12(s) = {D[0,1]}")
    print(f"D_21(s) = {D[1,0]}")
    print(f"D_22(s) = {D[1,1]}")

    print("\nMatrix N(s):")
    print(f"N_11(s) = {N[0,0]}")
    print(f"N_12(s) = {N[0,1]}")
    print(f"N_21(s) = {N[1,0]}")
    print(f"N_22(s) = {N[1,1]}")

solve_factorization()