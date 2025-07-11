import sympy as sp

def calculate_and_print_factorization():
    """
    This script prints the matrices for a left coprime factorization of H(s).
    """
    
    # Set up pretty printing for clear output
    sp.init_printing(use_unicode=True)
    
    # Define the symbolic variable s
    s = sp.symbols('s')

    # The left coprime factorization is H(s) = D(s)^-1 * N(s).
    # The derived matrices D(s) and N(s) are:
    D = sp.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sp.Matrix([
        [1, 1],
        [2, 0]
    ])

    print("The left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    
    print("\nD(s) =")
    sp.pprint(D)
    
    print("\nN(s) =")
    sp.pprint(N)
    
    print("\nIn equation form, this is:")
    print("\nH(s) = ")
    print(f"inv( [{D[0,0]}, {D[0,1]}] )   *   [{N[0,0]}, {N[0,1]}]")
    print(f"     [ {D[1,0]}, {D[1,1]}]           [ {N[1,0]}, {N[1,1]}]")

if __name__ == '__main__':
    calculate_and_print_factorization()