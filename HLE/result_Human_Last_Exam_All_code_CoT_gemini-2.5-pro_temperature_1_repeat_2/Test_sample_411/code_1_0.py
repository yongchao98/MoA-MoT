import sympy

def analyze_stc_diversity():
    """
    Analyzes the diversity order of three space-time codes using symbolic mathematics.
    """
    # --- Introduction to the method ---
    print("Analyzing the diversity order of the space-time codes.")
    print("The diversity order for a 2x1 MISO system is determined by the minimum rank of the codeword difference matrix, dS.")
    print("A code achieves full diversity (order 2) if det(dS) is non-zero for any non-zero error vector (e1, e2).")
    print("-" * 60)

    # --- Define symbolic variables ---
    # e1 and e2 represent the difference between symbols, e.g., e1 = x1 - x1'
    e1, e2 = sympy.symbols('e1 e2', complex=True)
    e1_c = sympy.conjugate(e1)
    e2_c = sympy.conjugate(e2)

    # --- Analysis for Code Sa ---
    print("(1) Analysis for Code Sa = [[x1, x2], [x2, x1]]")
    # Difference matrix for Sa
    dS_a = sympy.Matrix([[e1, e2], [e2, e1]])
    det_a = dS_a.det()
    print(f"The difference matrix dS_a is:\n{dS_a}")
    print(f"The determinant is: det(dS_a) = {det_a}")
    print("This determinant becomes zero if e1^2 - e2^2 = 0, which means e1 = e2 or e1 = -e2.")
    print("This condition can be met for non-zero e1 and e2. For such cases, the rank of dS_a is 1.")
    diversity_a = 1
    print(f"Conclusion: The minimum rank is 1. The diversity order for Sa is {diversity_a}.")
    print("-" * 60)

    # --- Analysis for Code Sb ---
    print("(2) Analysis for Code Sb = [[x1, x2], [x2, x1*]]")
    # Difference matrix for Sb
    dS_b = sympy.Matrix([[e1, e2], [e2, e1_c]])
    det_b = dS_b.det()
    # The symbolic output is e1*conjugate(e1) - e2**2, which is |e1|^2 - e2^2
    print(f"The difference matrix dS_b is:\n{dS_b}")
    print(f"The determinant is: det(dS_b) = {det_b}")
    print("This determinant becomes zero if |e1|^2 = e2^2. This is possible for non-zero e1 and e2.")
    print("For example, if e2 is a real number (e.g., e2=1) and e1 is a complex number with magnitude 1 (e.g., e1=1), the determinant is zero.")
    print("In this case, the rank of dS_b is 1.")
    diversity_b = 1
    print(f"Conclusion: The minimum rank is 1. The diversity order for Sb is {diversity_b}.")
    print("-" * 60)

    # --- Analysis for Code Sc ---
    print("(3) Analysis for Code Sc = [[-x1*, x2], [-x2*, -x1]]")
    # Difference matrix for Sc
    dS_c = sympy.Matrix([[-e1_c, e2], [-e2_c, -e1]])
    det_c = sympy.simplify(dS_c.det())
    # The symbolic output is e1*conjugate(e1) + e2*conjugate(e2), which is |e1|^2 + |e2|^2
    det_c_str = "|e1|^2 + |e2|^2"
    print(f"The difference matrix dS_c is:\n{dS_c}")
    print(f"The determinant is: det(dS_c) = {det_c} (i.e., {det_c_str})")
    print("This is a sum of non-negative terms. It is zero only if both |e1|^2=0 and |e2|^2=0, which implies e1=0 and e2=0.")
    print("Therefore, for any non-zero difference matrix dS_c, the determinant is non-zero, and its rank is 2.")
    diversity_c = 2
    print(f"Conclusion: The minimum rank is 2. The diversity order for Sc is {diversity_c}.")
    print("-" * 60)

    # --- Final Answers ---
    print("\n(a) Final Diversity Orders:")
    print(f"The diversity order for code Sa is: {diversity_a}")
    print(f"The diversity order for code Sb is: {diversity_b}")
    print(f"The diversity order for code Sc is: {diversity_c}")

    print("\n(b) Code with Maximum Diversity:")
    print("Assuming 'directivity' is a typo for 'diversity', code Sc provides the maximum diversity.")

if __name__ == '__main__':
    analyze_stc_diversity()