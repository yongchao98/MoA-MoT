import sympy as sp

def analyze_stbc():
    """
    Analyzes the diversity order and directivity of three space-time block codes.
    """
    # Define symbolic variables for the differences between symbols
    # d1 represents delta_x1 and d2 represents delta_x2
    d1 = sp.Symbol('delta_x1', complex=True, real=False)
    d2 = sp.Symbol('delta_x2', complex=True, real=False)

    # Define their conjugates
    d1_c = sp.conjugate(d1)
    d2_c = sp.conjugate(d2)

    print("Part (a): Diversity Order Analysis")
    print("="*40)
    print("The diversity order is the minimum rank of the codeword difference matrix, delta_S.")
    print("For a 2x2 matrix, the rank is 2 if the determinant is non-zero, and 1 if the determinant is zero (for a non-zero matrix).")
    print("We analyze the determinant of the difference matrix for each code.")
    print("The difference symbols delta_x1 and delta_x2 are not both zero for distinct codewords.\n")

    # --- Analysis for Code Sa ---
    print("--- Analysis for Code Sa ---")
    # The difference matrix is delta_S = [ [delta_x1, delta_x2], [delta_x2, delta_x1] ]
    delta_Sa = sp.Matrix([[d1, d2], [d2, d1]])
    det_delta_Sa = sp.det(delta_Sa)
    print("The difference matrix for Sa is:")
    sp.pprint(delta_Sa)
    print(f"\nThe determinant is: det(delta_Sa) = {det_delta_Sa}")
    print("This determinant becomes zero if delta_x1**2 - delta_x2**2 = 0. This can happen, for example, if we choose symbols such that delta_x1 = delta_x2.")
    print("Since the determinant can be zero for distinct codewords, the minimum rank is 1.")
    print("Final Equation: Diversity Order for Sa = 1\n")

    # --- Analysis for Code Sb ---
    print("--- Analysis for Code Sb ---")
    # The difference matrix is delta_S = [ [delta_x1, delta_x2], [delta_x2, conjugate(delta_x1)] ]
    delta_Sb = sp.Matrix([[d1, d2], [d2, d1_c]])
    det_delta_Sb = sp.det(delta_Sb)
    print("The difference matrix for Sb is:")
    sp.pprint(delta_Sb)
    print(f"\nThe determinant is: det(delta_Sb) = {det_delta_Sb}")
    print("This determinant becomes zero if Abs(delta_x1)**2 - delta_x2**2 = 0. This can happen if, for example, delta_x1 and delta_x2 are real-valued differences with equal magnitude.")
    print("Since the determinant can be zero, the minimum rank is 1.")
    print("Final Equation: Diversity Order for Sb = 1\n")

    # --- Analysis for Code Sc ---
    print("--- Analysis for Code Sc ---")
    # The difference matrix is delta_S = [ [-conjugate(delta_x1), delta_x2], [-conjugate(delta_x2), -delta_x1] ]
    delta_Sc = sp.Matrix([[-d1_c, d2], [-d2_c, -d1]])
    det_delta_Sc = sp.det(delta_Sc)
    print("The difference matrix for Sc is:")
    sp.pprint(delta_Sc)
    print(f"\nThe determinant is: det(delta_Sc) = {det_delta_Sc}")
    print("The term Abs(delta_x1)**2 + Abs(delta_x2)**2 is the sum of two non-negative numbers.")
    print("Since delta_x1 and delta_x2 are not both zero, this sum is always strictly positive.")
    print("The determinant is never zero, so the matrix is always full rank (rank 2).")
    print("Final Equation: Diversity Order for Sc = 2\n")

    # Part (b): Maximum Directivity
    print("Part (b): Maximum Directivity")
    print("="*40)
    print("In the context of space-time codes, 'directivity' is achieved through diversity. A higher diversity order provides a more reliable communication link.")
    print("Comparing the diversity orders:")
    print("  - Diversity(Sa) = 1")
    print("  - Diversity(Sb) = 1")
    print("  - Diversity(Sc) = 2")
    print("\nCode Sc achieves the maximum possible diversity order of 2 for a 2-antenna system. This is known as full diversity.")
    print("Therefore, Code Sc provides the maximum directivity.")

if __name__ == '__main__':
    analyze_stbc()