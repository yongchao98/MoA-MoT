import sympy

def solve_stc_analysis():
    """
    Analyzes the diversity order and performance of three given space-time codes.
    The script uses symbolic mathematics to derive the properties of each code
    and prints a detailed step-by-step explanation.
    """
    # Define symbolic variables for the error symbols e1 and e2.
    # e1 = x1 - x1', e2 = x2 - x2'. These are complex numbers.
    e1 = sympy.Symbol('e1', complex=True)
    e2 = sympy.Symbol('e2', complex=True)

    print("--- Analysis of Space-Time Codes ---")
    print("The diversity order of a space-time code is determined by the rank of the codeword difference matrix C.")
    print("For a 2x2 system, full diversity (order 2) is achieved if C is always full rank for any two distinct codewords.")
    print("This is equivalent to the determinant of C being non-zero for any non-zero error vector (e1, e2).\n")

    # --- Analysis for Part (a) ---
    print("Part (a): What is the diversity order for each code?\n")

    # --- Code Sa ---
    print("--- Analysis of Code Sa ---")
    Ca = sympy.Matrix([[e1, e2], [e2, e1]])
    det_Ca = Ca.det()
    print("The difference matrix Ca is:")
    sympy.pprint(Ca)
    print(f"\nThe determinant equation is: det(Ca) = {det_Ca} = 0")
    print("This can be written as: e1**2 - e2**2 = 0.")
    print("This determinant becomes zero if e1 = e2 or e1 = -e2. For example, if e1 = 2 and e2 = 2 (which are possible differences of QAM symbols),")
    print(f"the determinant is 2**2 - 2**2 = {2**2 - 2**2} = 0.")
    print("Since the determinant can be zero for non-zero error symbols, the matrix is not always full rank.")
    print("The minimum rank is 1. Therefore, the diversity order for Sa is 1.\n")

    # --- Code Sb ---
    print("--- Analysis of Code Sb ---")
    Cb = sympy.Matrix([[e1, e2], [e2, sympy.conjugate(e1)]])
    det_Cb = sympy.expand(Cb.det())
    print("The difference matrix Cb is:")
    sympy.pprint(Cb)
    # The result of det is e1*conjugate(e1) - e2**2, which is |e1|**2 - e2**2
    print(f"\nThe determinant equation is: det(Cb) = {det_Cb} = 0")
    print("This can be written as: |e1|**2 - e2**2 = 0.")
    print("This determinant can be zero. For this to happen, e2 must be a real number and |e1| must equal |e2|.")
    print("For example, if e1 = 2 and e2 = 2 (possible for QAM symbols),")
    print(f"the determinant is |2|**2 - 2**2 = {abs(2)**2 - 2**2} = 0.")
    print("Since the determinant can be zero, the matrix is not always full rank.")
    print("The minimum rank is 1. Therefore, the diversity order for Sb is 1.\n")

    # --- Code Sc ---
    print("--- Analysis of Code Sc ---")
    Cc = sympy.Matrix([[-sympy.conjugate(e1), e2], [-sympy.conjugate(e2), -e1]])
    det_Cc = sympy.expand(Cc.det())
    print("The difference matrix Cc is:")
    sympy.pprint(Cc)
    # The result is e1*conjugate(e1) + e2*conjugate(e2), which is |e1|^2 + |e2|^2
    print(f"\nThe determinant equation is: det(Cc) = {det_Cc} = 0")
    print("This can be written as: |e1|**2 + |e2|**2 = 0.")
    print("Since |e1|^2 and |e2|^2 are non-negative, their sum is zero only if both e1 = 0 and e2 = 0.")
    print("This corresponds to the case of identical codewords, which is excluded from the pairwise error analysis.")
    print("For any two distinct codewords, the determinant is non-zero, so the matrix is always full rank (rank 2).")
    print("Therefore, the diversity order for Sc is 2.\n")
    
    # --- Analysis for Part (b) ---
    print("Part (b): Which code provides the maximum directivity?\n")
    print("The term 'directivity' in this context is interpreted as the overall performance of the code.")
    print("Performance is primarily determined by the diversity order. A higher diversity order leads to significantly better error performance at high SNR.")
    print("Code Sc achieves the maximum possible diversity order of 2 (full diversity), while Sa and Sb only achieve a diversity order of 1.")
    print("Therefore, Sc provides the best performance or 'maximum directivity'.\n")

    print("<<<Diversity orders: Sa=1, Sb=1, Sc=2. Code providing maximum directivity: Sc>>>")

if __name__ == '__main__':
    solve_stc_analysis()