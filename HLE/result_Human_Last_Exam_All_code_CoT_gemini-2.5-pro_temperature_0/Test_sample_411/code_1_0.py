import sympy

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes (Sa, Sb, Sc)
    and determines which one provides the maximum diversity.
    """
    # Define symbolic variables for the complex differences of symbols.
    # delta_1 = x1 - e1, delta_2 = x2 - e2
    d1, d2 = sympy.symbols('delta_1 delta_2', complex=True)
    
    # Conjugates of the symbolic variables
    d1_c = sympy.conjugate(d1)
    d2_c = sympy.conjugate(d2)

    print("--- Analysis of Space-Time Block Codes ---")
    print("The diversity order is the minimum rank of the codeword difference matrix Delta_S.")
    print("For a 2x2 matrix, this means the determinant must be non-zero for all distinct codewords.\n")

    # --- (a) Diversity Order Analysis ---

    # Code Sa
    print("--- Analyzing Code Sa ---")
    Delta_Sa = sympy.Matrix([[d1, d2], [d2, d1]])
    det_Sa = sympy.simplify(Delta_Sa.det())
    print(f"Difference Matrix Delta_Sa:\n{Delta_Sa}")
    print(f"Determinant of Delta_Sa = {det_Sa}")
    print("Analysis: The determinant is zero if delta_1^2 = delta_2^2 (e.g., delta_1 = delta_2).")
    print("Since the determinant can be zero for non-zero differences, the matrix is not always full rank.")
    print("Diversity Order for Sa = 1\n")

    # Code Sb
    print("--- Analyzing Code Sb ---")
    Delta_Sb = sympy.Matrix([[d1, d2], [d2, d1_c]])
    det_Sb = sympy.simplify(Delta_Sb.det())
    print(f"Difference Matrix Delta_Sb:\n{Delta_Sb}")
    print(f"Determinant of Delta_Sb = {det_Sb}")
    print("Analysis: The determinant is |delta_1|^2 - delta_2^2. This can be zero for non-zero differences (e.g., delta_1=1, delta_2=1).")
    print("Since the determinant can be zero, the matrix is not always full rank.")
    print("Diversity Order for Sb = 1\n")

    # Code Sc
    print("--- Analyzing Code Sc ---")
    Delta_Sc = sympy.Matrix([[-d1_c, d2], [-d2_c, -d1]])
    det_Sc = sympy.simplify(Delta_Sc.det())
    print(f"Difference Matrix Delta_Sc:\n{Delta_Sc}")
    print(f"Determinant of Delta_Sc = {det_Sc}")
    print("Analysis: The determinant is |delta_1|^2 + |delta_2|^2.")
    print("Since at least one delta must be non-zero, this sum is always greater than 0.")
    print("The matrix is always full rank.")
    print("Diversity Order for Sc = 2 (Full Diversity)\n")

    # --- (b) Maximum Directivity ---
    print("--- (b) Maximum Directivity ---")
    print("Maximum 'directivity' in this context corresponds to the maximum diversity gain, which is achieved by the code with the highest diversity order.")
    print("Code Sc achieves the maximum possible diversity order of 2.\n")

    # --- Final Summary ---
    print("--- Summary of Results ---")
    print("(a) Diversity Orders:")
    print("    - Sa: 1")
    print("    - Sb: 1")
    print("    - Sc: 2")
    print("\n(b) Code with Maximum Directivity:")
    print("    - Sc")

if __name__ == '__main__':
    analyze_stbc()