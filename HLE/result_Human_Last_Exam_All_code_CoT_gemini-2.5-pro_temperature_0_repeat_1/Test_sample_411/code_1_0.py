import sympy as sp

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three space-time block codes (STBCs)
    for a 2x1 MISO system and determines which provides the maximum diversity.
    """
    # Define symbolic complex variables for the differences in transmitted symbols
    # dx1 = x1 - x1', dx2 = x2 - x2'
    dx1, dx2 = sp.symbols('dx1 dx2', complex=True)

    # --- Introduction ---
    print("--- Analysis of Space-Time Block Code Diversity Order ---")
    print("The diversity order of a space-time code is determined by the rank criterion.")
    print("A code has full diversity (order N=2 here) if the difference matrix dS = S - S'")
    print("has full rank for any two distinct codeword matrices S and S'.")
    print("For a 2x2 matrix dS, this is equivalent to checking if its determinant, det(dS),")
    print("is non-zero for any non-zero difference symbol vector (dx1, dx2).\n")

    # --- Code S_a ---
    print("--- (1) Analysis of Code S_a ---")
    # S_a = [[x1, x2], [x2, x1]]
    dS_a = sp.Matrix([[dx1, dx2], [dx2, dx1]])
    det_a = dS_a.det()
    print("The difference matrix for S_a is dS_a = S_a(x) - S_a(x'):")
    sp.pprint(dS_a)
    print("\nIts determinant is:")
    print(f"det(dS_a) = {det_a}")
    print("This determinant can be zero for non-zero dx1 and dx2.")
    print("For example, if dx1 = 1 and dx2 = 1, then det(dS_a) = 1**2 - 1**2 = 0.")
    print("Since the determinant can be zero, the matrix is not always full rank.")
    print("Result: The diversity order for S_a is 1.\n")

    # --- Code S_b ---
    print("--- (2) Analysis of Code S_b ---")
    # S_b = [[x1, x2], [x2, x1*]]
    dS_b = sp.Matrix([[dx1, dx2], [dx2, sp.conjugate(dx1)]])
    det_b = dS_b.det()
    print("The difference matrix for S_b is dS_b = S_b(x) - S_b(x'):")
    sp.pprint(dS_b)
    print("\nIts determinant is:")
    det_b_expr = sp.Abs(dx1)**2 - dx2**2
    print(f"det(dS_b) = |dx1|^2 - dx2^2")
    print("This determinant can also be zero for non-zero dx1 and dx2.")
    print("For example, if dx1 = 1 and dx2 = 1 (a real symbol), then det(dS_b) = |1|^2 - 1^2 = 0.")
    print("Since the determinant can be zero, the matrix is not always full rank.")
    print("Result: The diversity order for S_b is 1.\n")

    # --- Code S_c ---
    print("--- (3) Analysis of Code S_c ---")
    # S_c = [[-x1*, x2], [-x2*, -x1]]
    dS_c = sp.Matrix([[-sp.conjugate(dx1), dx2], [-sp.conjugate(dx2), -dx1]])
    det_c = dS_c.det()
    print("The difference matrix for S_c is dS_c = S_c(x) - S_c(x'):")
    sp.pprint(dS_c)
    print("\nIts determinant is:")
    det_c_expr = sp.Abs(dx1)**2 + sp.Abs(dx2)**2
    print(f"det(dS_c) = |dx1|^2 + |dx2|^2")
    print("This determinant is the sum of two non-negative terms, |dx1|^2 and |dx2|^2.")
    print("It is zero if and only if both dx1 = 0 and dx2 = 0.")
    print("Therefore, for any non-zero difference vector (dx1, dx2), the determinant is non-zero.")
    print("The matrix is always full rank, and the code achieves full diversity.")
    print("Result: The diversity order for S_c is 2.\n")

    # --- Summary and Conclusion ---
    print("--- Summary of Results ---")
    print("(a) What is the diversity order for each code?")
    print("    - Diversity order of S_a: 1")
    print("    - Diversity order of S_b: 1")
    print("    - Diversity order of S_c: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("Assuming 'directivity' refers to 'diversity', the code that provides the maximum")
    print("diversity is the one with the highest diversity order.")
    print("Therefore, S_c provides the maximum diversity.")

if __name__ == '__main__':
    analyze_stbc_diversity()