import sympy
from sympy import symbols, Matrix, conjugate, det, pprint

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes.
    """
    # Define symbolic variables for two distinct codewords
    x1, x2, x1_p, x2_p = symbols('x1 x2 x1_p x2_p', complex=True)

    # Define the error symbols e1 and e2
    e1 = x1 - x1_p
    e2 = x2 - x2_p
    
    # --- Analysis for Code Sa ---
    print("--- Analysis for Code Sa ---")
    print("S_a = [x1, x2]")
    print("      [x2, x1]\n")
    
    # Define the difference matrix Delta_Sa
    # Delta_Sa = S_a(x1, x2) - S_a(x1_p, x2_p)
    Delta_Sa = Matrix([[e1, e2], [e2, e1]])
    
    print("The codeword difference matrix is Delta_S_a = S(x) - S(x'):")
    pprint(Delta_Sa)
    print()

    # Calculate the determinant
    det_Sa = det(Delta_Sa)
    
    print("The determinant of the difference matrix is:")
    print(f"det(Delta_S_a) = {det_Sa}")
    print("\nAnalysis:")
    print("The determinant is e1**2 - e2**2. This expression can be zero for non-zero error symbols.")
    print("For example, if e1 = e2 (e.g., x1-x1' = x2-x2' != 0), the determinant is 0.")
    print("Since the difference matrix can be rank-deficient for some error events, the code does not provide full diversity.")
    print("The diversity order for Code Sa is 1.")
    print("-" * 40)

    # --- Analysis for Code Sb ---
    print("--- Analysis for Code Sb ---")
    print("S_b = [x1,   x2]")
    print("      [x2, x1*]\n")

    # Define the difference matrix Delta_Sb
    # Delta_Sb = S_b(x1, x2) - S_b(x1_p, x2_p)
    # The difference of conjugates is the conjugate of the difference: x1* - x1_p* = (x1 - x1_p)* = e1*
    Delta_Sb = Matrix([[e1, e2], [e2, conjugate(e1)]])

    print("The codeword difference matrix is Delta_S_b = S(x) - S(x'):")
    pprint(Delta_Sb)
    print()

    # Calculate the determinant
    det_Sb = det(Delta_Sb)

    print("The determinant of the difference matrix is:")
    print(f"det(Delta_S_b) = {det_Sb}")
    print("This can be written as |e1|^2 - e2^2.")
    print("\nAnalysis:")
    print("This determinant can also be zero for non-zero error symbols.")
    print("For example, if the symbols are from a real constellation (like BPSK), choosing e1=e2 makes the determinant zero.")
    print("For a complex QAM constellation, we can also find cases where |e1|^2 = e2^2 (e.g., e1 = 3+4j, e2 = 5).")
    print("Since the difference matrix can be rank-deficient, the code does not provide full diversity.")
    print("The diversity order for Code Sb is 1.")
    print("-" * 40)
    
    # --- Analysis for Code Sc ---
    print("--- Analysis for Code Sc ---")
    print("S_c = [-x1*,  x2]")
    print("      [-x2*, -x1]\n")
    
    # Define the difference matrix Delta_Sc
    # Delta_Sc = S_c(x1, x2) - S_c(x1_p, x2_p)
    Delta_Sc = Matrix([[-conjugate(e1), e2], [-conjugate(e2), -e1]])

    print("The codeword difference matrix is Delta_S_c = S(x) - S(x'):")
    pprint(Delta_Sc)
    print()
    
    # Calculate the determinant
    det_Sc = det(Delta_Sc).simplify()

    print("The determinant of the difference matrix is:")
    # The result of det() is e1*conjugate(e1) + e2*conjugate(e2)
    # Let's show it in a more familiar form |e1|^2 + |e2|^2
    e1_abs_sq = symbols('|e1|^2')
    e2_abs_sq = symbols('|e2|^2')
    print(f"det(Delta_S_c) = {det_Sc} = {e1_abs_sq} + {e2_abs_sq}")
    print("\nAnalysis:")
    print("The determinant is the sum of the squared magnitudes of the error symbols.")
    print("This sum is zero ONLY if both e1 and e2 are zero, which means there is no error (x=x').")
    print("For any non-zero error vector e (i.e., at least one of e1, e2 is non-zero), the determinant is strictly positive.")
    print("Therefore, the difference matrix is always full-rank, and the code achieves full diversity.")
    print("The diversity order for Code Sc is 2.")
    print("-" * 40)

    # --- Final Conclusion ---
    print("\n(a) Summary of Diversity Orders:")
    print("Diversity order of Sa: 1")
    print("Diversity order of Sb: 1")
    print("Diversity order of Sc: 2")
    
    print("\n(b) Code with Maximum Diversity:")
    print("Code Sc provides the maximum diversity order (which is the full diversity order of 2 for this system).")


if __name__ == '__main__':
    analyze_stbc()