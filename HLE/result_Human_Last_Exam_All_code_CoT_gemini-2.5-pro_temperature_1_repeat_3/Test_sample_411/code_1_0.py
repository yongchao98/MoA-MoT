import sympy

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three 2x2 space-time block codes (STBCs)
    and determines which one provides the maximum diversity.
    """
    print("This script analyzes the diversity order of three 2x2 space-time block codes.")
    print("The diversity order is determined by the Rank-Diversity Criterion.")
    print("This criterion states that the diversity order is the minimum rank of the codeword difference")
    print("matrix Delta(e) = S(x) - S(x') for any non-zero error vector e = x - x'.")
    print("For a 2x2 code, full diversity (order 2) is achieved if the determinant of the")
    print("difference matrix is non-zero for any non-zero error vector.\n")

    # Define symbolic complex variables for the error vector e = [e1, e2]
    e1 = sympy.Symbol('e1', complex=True)
    e2 = sympy.Symbol('e2', complex=True)

    # --- Analysis of Code S_a ---
    print("--- Analyzing Code S_a = [[x1, x2], [x2, x1]] ---")
    delta_a = sympy.Matrix([[e1, e2], [e2, e1]])
    print("The difference matrix Delta_a is:")
    sympy.pprint(delta_a)
    
    det_a = delta_a.det()
    print(f"\nThe determinant of Delta_a is: det(Delta_a) = {det_a}")
    
    print("\nWe check if the determinant can be zero for a non-zero error vector (e1, e2) != (0, 0).")
    print("det(Delta_a) = 0  =>  e1**2 = e2**2  =>  e1 = e2 or e1 = -e2.")
    print("For any M-QAM constellation, it's possible to find distinct symbols x1, x'1, x2, x'2")
    print("such that the differences e1 = x1 - x'1 and e2 = x2 - x'2 satisfy e1 = e2 != 0.")
    print("For instance, in a 16-QAM system, we can have e1 = 2 and e2 = 2.")
    print("For this non-zero error vector, det(Delta_a) = 2**2 - 2**2 = 0.")
    print("Since the determinant can be zero, the matrix can be rank-deficient (rank 1).")
    diversity_a = 1
    print(f"Therefore, the minimum rank is 1, and the diversity order for code S_a is {diversity_a}.\n")

    # --- Analysis of Code S_b ---
    print("--- Analyzing Code S_b = [[x1, x2], [x2, x1*]] ---")
    delta_b = sympy.Matrix([[e1, e2], [e2, sympy.conjugate(e1)]])
    print("The difference matrix Delta_b is:")
    sympy.pprint(delta_b)

    det_b = delta_b.det()
    print(f"\nThe determinant of Delta_b is: det(Delta_b) = {det_b}")
    
    print("\nWe check if the determinant can be zero for a non-zero error vector (e1, e2) != (0, 0).")
    print("det(Delta_b) = 0  =>  e1*conjugate(e1) = e2**2  =>  |e1|**2 = e2**2.")
    print("For this equality to hold with a complex e2, the imaginary part of e2**2 must be zero,")
    print("which means e2 must be a real or purely imaginary number.")
    print("Let's assume e2 is real. We check if we can find e1 (complex) and e2 (real) in the")
    print("QAM difference set such that |e1| = |e2|.")
    print("Using a 16-QAM difference set, we can find e2 = 4 (real) and e1 = 4j (complex).")
    print(f"For this non-zero error vector, det(Delta_b) = |4j|^2 - 4^2 = {abs(4j)**2} - {4**2} = 16 - 16 = 0.")
    print("Since the determinant can be zero, the matrix can be rank-deficient (rank 1).")
    diversity_b = 1
    print(f"Therefore, the minimum rank is 1, and the diversity order for code S_b is {diversity_b}.\n")

    # --- Analysis of Code S_c ---
    print("--- Analyzing Code S_c = [[-x1*, x2], [-x2*, -x1]] ---")
    delta_c = sympy.Matrix([[-sympy.conjugate(e1), e2], [-sympy.conjugate(e2), -e1]])
    print("The difference matrix Delta_c is:")
    sympy.pprint(delta_c)

    det_c = delta_c.det()
    # sympy.expand is used for a cleaner output
    print(f"\nThe determinant of Delta_c is: det(Delta_c) = {sympy.expand(det_c)}")
    
    print("\nWe check if the determinant can be zero for a non-zero error vector (e1, e2) != (0, 0).")
    print("det(Delta_c) = 0  =>  e1*conjugate(e1) + e2*conjugate(e2) = 0  =>  |e1|**2 + |e2|**2 = 0.")
    print("Since |e1|**2 >= 0 and |e2|**2 >= 0, their sum is zero if and only if |e1|**2 = 0 and |e2|**2 = 0.")
    print("This implies e1 = 0 and e2 = 0, which corresponds to the zero error vector.")
    print("For any *non-zero* error vector, the determinant is strictly positive and thus never zero.")
    print("The matrix Delta_c is always full rank (rank 2) for any non-zero error vector.")
    diversity_c = 2
    print(f"Therefore, the minimum rank is 2, and the diversity order for code S_c is {diversity_c}.\n")
    
    # --- Final Conclusion ---
    print("--- Conclusion ---")
    print("(a) What is the diversity order for each code?")
    print(f"The diversity order for S_a is {diversity_a}.")
    print(f"The diversity order for S_b is {diversity_b}.")
    print(f"The diversity order for S_c is {diversity_c}.")

    print("\n(b) Which code provides the maximum directivity (diversity)?")
    print(f"Code S_c provides the maximum diversity order of {diversity_c}, which is the full diversity for a 2-antenna system.")

if __name__ == '__main__':
    analyze_stbc_diversity()