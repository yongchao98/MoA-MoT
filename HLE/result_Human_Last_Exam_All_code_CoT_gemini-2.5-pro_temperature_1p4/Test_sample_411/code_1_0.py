import sympy

def analyze_space_time_codes():
    """
    Analyzes three space-time codes to determine their diversity order.
    """
    # Define symbolic variables for two distinct symbol vectors x=(x1, x2) and e=(e1, e2)
    x1, x2, e1, e2 = sympy.symbols('x1 x2 e1 e2', complex=True)
    
    # Define the difference symbols
    d1 = x1 - e1
    d2 = x2 - e2

    print("Diversity analysis for the given space-time codes.\n")
    print("Let d1 = x1 - e1 and d2 = x2 - e2.")
    print("A code has full diversity (order 2) if det(Delta_S) is non-zero whenever (d1, d2) != (0, 0).\n")
    print("-" * 60)

    # --- Analysis for Code Sa ---
    print("Part (a): Diversity Order Analysis\n")
    print("1. Analysis for Code Sa")
    S_a = sympy.Matrix([[x1, x2], [x2, x1]])
    Delta_S_a = sympy.Matrix([[d1, d2], [d2, d1]])
    det_S_a = Delta_S_a.det()
    
    print("S_a = ")
    sympy.pprint(S_a)
    print("\nDelta_S_a = S(x1, x2) - S(e1, e2) =")
    sympy.pprint(Delta_S_a)
    print(f"\ndet(Delta_S_a) = {det_S_a}")
    print("This determinant can be zero if d1^2 = d2^2 (e.g., if d1 = d2), even if d1 and d2 are not zero.")
    print("Therefore, the code is not full diversity.")
    print("Diversity Order for S_a = 1\n")
    print("-" * 60)

    # --- Analysis for Code Sb ---
    print("2. Analysis for Code Sb")
    S_b = sympy.Matrix([[x1, x2], [x2, sympy.conjugate(x1)]])
    Delta_S_b = S_b.subs({x1:d1, x2:d2}) # This is S_b(d1,d2)
    Delta_S_b_true = sympy.Matrix([[d1, d2], [d2, sympy.conjugate(d1)]])
    det_S_b = Delta_S_b_true.det()
    
    print("S_b = ")
    sympy.pprint(S_b)
    print("\nDelta_S_b = S(x1, x2) - S(e1, e2) =")
    sympy.pprint(Delta_S_b_true)
    print(f"\ndet(Delta_S_b) = {sympy.simplify(det_S_b)}")
    print("This determinant can be zero if |d1|^2 = d2^2. This is possible for non-zero d1, d2 if d2 is real and |d1| = |d2|.")
    print("Therefore, the code is not full diversity.")
    print("Diversity Order for S_b = 1\n")
    print("-" * 60)
    
    # --- Analysis for Code Sc ---
    print("3. Analysis for Code Sc")
    S_c = sympy.Matrix([[-sympy.conjugate(x1), x2], [-sympy.conjugate(x2), -x1]])
    Delta_S_c = sympy.Matrix([[-sympy.conjugate(d1), d2], [-sympy.conjugate(d2), -d1]])
    det_S_c = Delta_S_c.det()
    
    print("S_c = ")
    sympy.pprint(S_c)
    print("\nDelta_S_c = S(x1, x2) - S(e1, e2) =")
    sympy.pprint(Delta_S_c)
    print(f"\ndet(Delta_S_c) = {sympy.simplify(det_S_c)}")
    print("Since d1*conjugate(d1) = |d1|^2 and d2*conjugate(d2) = |d2|^2, the determinant is |d1|^2 + |d2|^2.")
    print("This is zero only if d1 = 0 and d2 = 0, which means the codewords are identical.")
    print("The determinant is always non-zero for distinct codewords.")
    print("Therefore, the code is full diversity.")
    print("Diversity Order for S_c = 2\n")
    print("-" * 60)

    # --- Conclusion for Part (b) ---
    print("Part (b): Code with Maximum Directivity\n")
    print("Interpreting 'directivity' as 'diversity', the code with the maximum diversity order is S_c.")
    print("It achieves the full diversity of 2 for a 2x1 MISO system.")

if __name__ == '__main__':
    analyze_space_time_codes()
