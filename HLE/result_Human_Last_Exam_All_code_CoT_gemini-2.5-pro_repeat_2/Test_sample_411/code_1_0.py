import sympy

def analyze_stbc():
    """
    Analyzes the diversity order of three 2x2 STBCs.
    """
    # Define symbolic variables for the differences of complex symbols
    # delta_1 = x_1 - x'_1, delta_2 = x_2 - x'_2
    d1, d2 = sympy.symbols('delta_1, delta_2', complex=True)

    print("The diversity order of an STBC is determined by the rank of the difference matrix Delta_S.")
    print("For full diversity of 2, the determinant of Delta_S must be non-zero for any two distinct codewords.")
    print("Let the difference symbols be delta_1 and delta_2, which cannot both be zero.\n")

    # --- Analysis for Code S_a ---
    print("--- Analysis for Code S_a = [[x1, x2], [x2, x1]] ---")
    # Difference matrix for Sa
    Delta_Sa = sympy.Matrix([[d1, d2], [d2, d1]])
    # Determinant of the difference matrix
    det_Sa = Delta_Sa.det()
    
    print("The difference matrix Delta_S_a is:")
    print(sympy.pretty(Delta_Sa, use_unicode=False))
    print("\nIts determinant is:")
    # Print the equation for the determinant
    print(f"det(Delta_S_a) = ({d1})*({d1}) - ({d2})*({d2})")
    print(f"              = {det_Sa}")
    print("\nThis determinant equals zero if delta_1^2 = delta_2^2 (e.g., if delta_1 = delta_2).")
    print("Since this can happen for non-zero delta_1 and delta_2, the matrix can be rank-deficient.")
    print("Result: The diversity order for code S_a is 1.")
    print("-" * 50)

    # --- Analysis for Code S_b ---
    print("\n--- Analysis for Code S_b = [[x1, x2], [x2, x1*]] ---")
    # Difference matrix for Sb
    Delta_Sb_str = f"[[{d1}, {d2}], [{d2}, conjugate({d1})]]"
    
    print("The difference matrix Delta_S_b is:")
    print(Delta_Sb_str)
    print("\nIts determinant is:")
    # Print the equation for the determinant
    print(f"det(Delta_S_b) = ({d1})*(conjugate({d1})) - ({d2})*({d2})")
    print(f"              = |{d1}|^2 - ({d2})^2")
    print("\nThis determinant can be zero for non-zero delta_1 and delta_2.")
    print("For example, if the symbols are from a BPSK constellation, delta_1 and delta_2 will be real.")
    print("If delta_1 = 2 and delta_2 = 2, the determinant is |2|^2 - 2^2 = 0.")
    print("Thus, the matrix can be rank-deficient.")
    print("Result: The diversity order for code S_b is 1.")
    print("-" * 50)

    # --- Analysis for Code S_c ---
    print("\n--- Analysis for Code S_c = [[-x1*, x2], [-x2*, -x1]] ---")
    # Difference matrix for Sc
    Delta_Sc_str = f"[[-conjugate({d1}), {d2}], [-conjugate({d2}), -{d1}]]"
    
    print("The difference matrix Delta_S_c is:")
    print(Delta_Sc_str)
    print("\nIts determinant is:")
    # Print the equation for the determinant
    print(f"det(Delta_S_c) = (-conjugate({d1}))*(-{d1}) - ({d2})*(-conjugate({d2}))")
    print(f"              = |{d1}|^2 + |{d2}|^2")
    print("\nThis determinant is the sum of squares of magnitudes.")
    print("It is zero only if both delta_1 = 0 and delta_2 = 0, which implies the codewords are identical.")
    print("For any two distinct codewords, the determinant is always non-zero.")
    print("Thus, the matrix is always full-rank.")
    print("Result: The diversity order for code S_c is 2.")
    print("-" * 50)

    print("\n--- Summary ---")
    print("(a) The diversity orders are:")
    print("    - Code S_a: 1")
    print("    - Code S_b: 1")
    print("    - Code S_c: 2")
    print("\n(b) Which code provides the maximum diversity?")
    print("    The maximum possible diversity for this system is 2.")
    print("    Code S_c is the only one that achieves this maximum diversity.")

if __name__ == '__main__':
    analyze_stbc()