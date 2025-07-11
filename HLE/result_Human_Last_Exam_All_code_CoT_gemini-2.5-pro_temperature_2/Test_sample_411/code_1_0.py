import sympy
from sympy import symbols, Matrix, conjugate, I

def analyze_diversity():
    """
    Analyzes the diversity order of three space-time codes S_a, S_b, and S_c.
    """
    # Define complex symbols for the error vector components delta_1 and delta_2
    # delta_1 = x1 - e1, delta_2 = x2 - e2
    d1, d2 = symbols('delta_1, delta_2', complex=True)

    # --- Introduction ---
    print("Analyzing diversity order using the Rank-Determinant Criterion.")
    print("A code has full diversity order (2) if the determinant of its difference matrix is non-zero for any two distinct symbol vectors.")
    print("The difference symbols are delta_1 = x1 - e1 and delta_2 = x2 - e2.\n")

    # --- Part (a): Diversity Order Analysis ---
    print("Part (a): Diversity Order for each code")
    print("="*40)

    # --- Analysis for Code S_a ---
    print("\nAnalysis for Code S_a = [[x1, x2], [x2, x1]]")
    print("-" * 40)
    delta_Sa = Matrix([[d1, d2], [d2, d1]])
    det_Sa = sympy.simplify(delta_Sa.det())
    print(f"The difference matrix is:\n{delta_Sa}\n")
    print(f"The determinant is: det(ΔS_a) = {det_Sa}")
    print("\nAnalysis:")
    print("The determinant is zero if delta_1**2 = delta_2**2, which means delta_1 = delta_2 or delta_1 = -delta_2.")
    print("This condition can be met for non-zero error symbols. For example, using BPSK symbols {1, -1}:")
    print("if x = (1, 1) and e = (-1, -1), then delta_1 = 2 and delta_2 = 2, making the determinant 0.")
    print("Since the determinant can be zero, the difference matrix can be rank-deficient (rank 1).")
    print("Diversity Order for S_a = 1")

    # --- Analysis for Code S_b ---
    print("\nAnalysis for Code S_b = [[x1, x2], [x2, x1*]]")
    print("-" * 40)
    delta_Sb = Matrix([[d1, d2], [d2, conjugate(d1)]])
    det_Sb = sympy.simplify(delta_Sb.det())
    print(f"The difference matrix is:\n{delta_Sb}\n")
    # Using specific terms for clarity in the final equation output
    print(f"The determinant is: det(ΔS_b) = delta_1 * conjugate(delta_1) - delta_2**2 = |delta_1|**2 - delta_2**2")
    print("\nAnalysis:")
    print("The determinant is zero if |delta_1|**2 = delta_2**2.")
    print("This requires the imaginary part of delta_2**2 to be zero, meaning delta_2 must be purely real or purely imaginary.")
    print("For constellations like 64-QAM, we can find symbol differences that satisfy this.")
    print("For example, based on the Pythagorean triple (6, 8, 10), we can find symbol differences for 64-QAM such that delta_1 = 6 + 8j and delta_2 = 10.")
    print("In this case, |delta_1|**2 = 6**2 + 8**2 = 100, and delta_2**2 = 10**2 = 100. The determinant is zero.")
    print("Since the determinant can be zero, the difference matrix can be rank-deficient (rank 1).")
    print("Diversity Order for S_b = 1")


    # --- Analysis for Code S_c ---
    print("\nAnalysis for Code S_c = [[-x1*, x2], [-x2*, -x1]]")
    print("-" * 40)
    delta_Sc = Matrix([[-conjugate(d1), d2], [-conjugate(d2), -d1]])
    det_Sc = sympy.simplify(delta_Sc.det())
    print(f"The difference matrix is:\n{delta_Sc}\n")
    print(f"The determinant is: det(ΔS_c) = (-conjugate(delta_1))*(-delta_1) - (d2)*(-conjugate(d2)) = |delta_1|**2 + |delta_2|**2")
    print("\nAnalysis:")
    print("The determinant is the sum of two non-negative terms, |delta_1|**2 and |delta_2|**2.")
    print("This sum is zero if and only if both terms are zero, i.e., |delta_1| = 0 and |delta_2| = 0.")
    print("This implies delta_1 = 0 and delta_2 = 0, which corresponds to the case where the transmitted symbol vectors are identical (x=e), not distinct.")
    print("For any two distinct codewords, the determinant is always non-zero. The matrix is always full rank (2).")
    print("Diversity Order for S_c = 2")

    # --- Part (b): Conclusion on Maximum Diversity ---
    print("\n" + "="*40)
    print("Part (b): Which code provides the maximum diversity?")
    print("="*40)
    print("\nThe diversity orders are 1, 1, and 2 for S_a, S_b, and S_c, respectively.")
    print("The maximum possible diversity for a 2x1 system is 2.")
    print("Code S_c achieves this maximum diversity order.")
    print("\nConclusion: Code S_c provides the maximum diversity.")

if __name__ == "__main__":
    analyze_diversity()