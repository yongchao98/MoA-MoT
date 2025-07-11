import sympy as sp

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes Sa, Sb, and Sc.
    """
    # Define e1 and e2 as complex symbols for the error vector [e1, e2]
    e1, e2 = sp.symbols('e1 e2', complex=True)

    print("This script analyzes the diversity order of three 2x2 STBCs.")
    print("The diversity order is the minimum rank of the codeword difference matrix.")
    print("For a 2x2 matrix, this means checking if its determinant can be zero for non-zero errors.\n")

    # --- Analysis for Code Sa ---
    print("--- (1) Analysis for Code Sa ---")
    # S_a = [[x1, x2], [x2, x1]]
    delta_Sa = sp.Matrix([[e1, e2], [e2, e1]])
    det_Sa = sp.det(delta_Sa)
    
    print("The difference matrix Delta_S_a is:")
    sp.pprint(delta_Sa)
    print(f"\nThe determinant is: {det_Sa}")
    print("This determinant, e1**2 - e2**2, can be zero for non-zero errors.")
    print("For example, if e1 = 1 and e2 = 1 (a valid non-zero error), det(Delta_S_a) = 1**2 - 1**2 = 0.")
    print("When the determinant is zero, the matrix rank is less than 2. In this case, it is 1.")
    print("The diversity order for Sa is 1.")
    print("\n" + "="*40 + "\n")

    # --- Analysis for Code Sb ---
    print("--- (2) Analysis for Code Sb ---")
    # S_b = [[x1, x2], [x2, x1*]]
    delta_Sb = sp.Matrix([[e1, e2], [e2, sp.conjugate(e1)]])
    det_Sb = sp.det(delta_Sb)

    print("The difference matrix Delta_S_b is:")
    sp.pprint(delta_Sb)
    print(f"\nThe determinant is: {det_Sb}")
    print("This determinant, e1*conjugate(e1) - e2**2 or |e1|**2 - e2**2, can also be zero for non-zero errors.")
    print("For example, if e1 = 1 and e2 = 1, det(Delta_S_b) = |1|**2 - 1**2 = 0.")
    print("When the determinant is zero, the matrix rank is 1.")
    print("The diversity order for Sb is 1.")
    print("\n" + "="*40 + "\n")

    # --- Analysis for Code Sc ---
    print("--- (3) Analysis for Code Sc ---")
    # S_c = [[-x1*, x2], [-x2*, -x1]]
    delta_Sc = sp.Matrix([[-sp.conjugate(e1), e2], [-sp.conjugate(e2), -e1]])
    det_Sc = sp.det(delta_Sc)

    print("The difference matrix Delta_S_c is:")
    sp.pprint(delta_Sc)
    print(f"\nThe determinant is: {sp.simplify(det_Sc)}")
    print("This determinant, e1*conjugate(e1) + e2*conjugate(e2) or |e1|**2 + |e2|**2, is the sum of two non-negative terms.")
    print("This expression is zero if and only if both |e1|**2 = 0 and |e2|**2 = 0, which means e1 = 0 and e2 = 0.")
    print("For any non-zero error vector, the determinant is non-zero, so the matrix is always full rank (rank=2).")
    print("The diversity order for Sc is 2.")
    print("\n" + "="*40 + "\n")
    
    # --- Final Conclusion ---
    print("--- (a) Conclusion on Diversity Order ---")
    print("The diversity order for code Sa is 1.")
    print("The diversity order for code Sb is 1.")
    print("The diversity order for code Sc is 2.")
    
    print("\n--- (b) Conclusion on Maximum Directivity ---")
    print("Maximum 'directivity' corresponds to the maximum diversity gain, achieved by the code with the highest diversity order.")
    print("Code Sc provides the maximum directivity because its diversity order is 2, which is the maximum possible for a 2x1 system.")

if __name__ == '__main__':
    # To run this code, you need to have sympy installed (`pip install sympy`)
    analyze_stbc()