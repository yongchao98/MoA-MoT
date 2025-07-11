import sympy

def main():
    """
    Calculates and displays a left coprime factorization of the given transfer function.
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')
    sympy.init_printing(use_unicode=True)

    # The left coprime factorization is H(s) = D(s)^-1 * N(s).
    # Based on the step-by-step derivation, the polynomial matrices D(s) and N(s) are:

    # The denominator matrix D(s)
    D = sympy.Matrix([[1, s - 1], [0, s**2 - 1]])

    # The numerator matrix N(s)
    N = sympy.Matrix([[1, 1], [2, 0]])

    print("A left coprime factorization of H(s) is given by H(s) = D(s)^-1 * N(s), where:")
    
    print("\nD(s) =")
    sympy.pprint(D)

    print("\nN(s) =")
    sympy.pprint(N)
    
    # Display the full equation H(s) = D(s)^-1 * N(s) with each part shown
    print("\nThe full factorization equation is:")
    
    # Calculate D^-1(s) to display
    D_inv = sympy.simplify(D.inv())
    
    # Calculate H(s) = D^-1(s) * N(s)
    H_calc = sympy.simplify(D_inv * N)

    print("\nH(s) = ")
    sympy.pprint(H_calc)
    
    print("\nwhich is equal to ( D(s)^-1 ) * N(s) :")
    
    print("\nD(s)^-1 =")
    sympy.pprint(D_inv)

    print("\n*")

    print("\nN(s) =")
    sympy.pprint(N)

if __name__ == "__main__":
    main()
