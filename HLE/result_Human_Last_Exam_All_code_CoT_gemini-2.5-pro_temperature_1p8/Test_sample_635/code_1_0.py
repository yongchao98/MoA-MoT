import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and verifies a left coprime factorization for the given transfer function H(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # Based on our derivation, we propose the following matrices D(s) and N(s)
    D = sympy.Matrix([[1, s - 1],
                      [0, s**2 - 1]])

    N = sympy.Matrix([[1, 1],
                      [2, 0]])

    # The original transfer function H(s)
    H = sympy.Matrix([[(s - 1) / (s + 1), 1],
                      [2 / (s**2 - 1), 0]])

    print("Calculating a left coprime factorization H(s) = D(s)^-1 * N(s) ...\n")
    
    # --- Print the final answer in the required format ---
    print("The resulting matrices are:")
    print("D(s) =")
    print(f"[[ {D[0,0]} , {D[0,1]} ],")
    print(f" [ {D[1,0]} , {D[1,1]} ]]")
    
    print("\nN(s) =")
    print(f"[[ {N[0,0]} , {N[0,1]} ],")
    print(f" [ {N[1,0]} , {N[1,1]} ]]")
    print("\n" + "="*40 + "\n")

    # --- Verification Steps ---
    print("Verification steps:")
    
    # 1. Verify the factorization D(s)^-1 * N(s) = H(s)
    print("1. Checking if D(s)^-1 * N(s) equals H(s)...")
    try:
        D_inv = D.inv()
        H_calculated = sympy.simplify(D_inv * N)
        
        print("Original H(s):")
        sympy.pretty_print(H)
        
        print("\nCalculated D(s)^-1 * N(s):")
        sympy.pretty_print(H_calculated)

        if H_calculated == H:
            print("\nResult: The factorization is correct.")
        else:
            print("\nResult: The factorization is incorrect.")
    except Exception as e:
        print(f"An error occurred during verification: {e}")

    # 2. Verify coprimality by checking the rank of [D(s) N(s)]
    print("\n2. Checking if D(s) and N(s) are left coprime...")
    M = D.row_join(N)
    print("The composite matrix [D(s) N(s)] is:")
    sympy.pretty_print(M)
    
    # Find a 2x2 minor with a non-zero constant determinant
    minor_13 = M[:, [0, 2]] # Columns 1 and 3 of M
    det_minor_13 = minor_13.det()
    
    print(f"\nThe determinant of the minor formed by columns 1 and 3 is: {det_minor_13}")
    
    if det_minor_13 != 0 and det_minor_13.is_constant():
        print("Result: Since the determinant is a non-zero constant, [D(s) N(s)] has full rank for all s.")
        print("Therefore, D(s) and N(s) are left coprime.")
    else:
        print("Result: The matrices may not be left coprime as no simple constant minor was found.")

if __name__ == '__main__':
    calculate_left_coprime_factorization()