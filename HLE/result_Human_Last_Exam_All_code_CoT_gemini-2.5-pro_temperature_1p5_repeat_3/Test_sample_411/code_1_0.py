import sympy

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three space-time block codes.
    """
    # Define complex symbol differences, delta_1 = x1 - x'1 and delta_2 = x2 - x'2.
    # For distinct codewords, at least one of delta_1 or delta_2 must be non-zero.
    d1 = sympy.Symbol('delta_1', complex=True)
    d2 = sympy.Symbol('delta_2', complex=True)

    # --- Part (a): Find the diversity order for each code ---
    print("(a) What is the diversity order for each code?")
    print("The diversity order depends on the rank of the codeword difference matrix, Delta_S.")
    print("If det(Delta_S) is non-zero for all distinct codewords, the diversity order is 2.")
    print("If det(Delta_S) can be zero for some distinct codewords, the diversity order is 1.")

    # --- Analysis for Code Sa ---
    # S_a = [[x1, x2], [x2, x1]]
    # Delta_S_a = [[d1, d2], [d2, d1]]
    det_Sa = d1**2 - d2**2
    print("\n--- Code Sa ---")
    print(f"The determinant of the difference matrix is: det(Delta_S_a) = delta_1**2 - delta_2**2")
    print("This determinant is zero if delta_1 = delta_2 or delta_1 = -delta_2.")
    print("This condition can be met for distinct codewords (non-zero delta_1, delta_2).")
    print("Therefore, the minimum rank is 1, and the diversity order is 1.")

    # --- Analysis for Code Sb ---
    # S_b = [[x1, x2], [x2, x1*]]
    # Delta_S_b = [[d1, d2], [d2, d1.conjugate()]]
    det_Sb_expr = "|delta_1|^2 - delta_2**2"
    print("\n--- Code Sb ---")
    print(f"The determinant of the difference matrix is: det(Delta_S_b) = {det_Sb_expr}")
    print("This determinant is zero if |delta_1|^2 = delta_2**2.")
    print("This condition can also be met for distinct codewords (e.g., using QAM symbols).")
    print("Therefore, the minimum rank is 1, and the diversity order is 1.")

    # --- Analysis for Code Sc ---
    # S_c = [[-x1*, x2], [-x2*, -x1]]
    # Delta_S_c = [[-d1*, d2], [-d2*, -d1]]
    det_Sc_expr = "|delta_1|^2 + |delta_2|^2"
    print("\n--- Code Sc ---")
    print(f"The determinant of the difference matrix is: det(Delta_S_c) = {det_Sc_expr}")
    print("This determinant is a sum of non-negative terms. It is zero only if both delta_1 = 0 and delta_2 = 0.")
    print("This only happens when the codewords are identical, not distinct.")
    print("Therefore, the determinant is always non-zero for distinct codewords, the matrix is always full-rank, and the diversity order is 2.")

    # --- Part (b): Determine which code has maximum diversity ---
    print("\n(b) Which code provides the maximum directivity?")
    print("Assuming 'directivity' refers to 'diversity', Code Sc provides the maximum diversity order.")

if __name__ == '__main__':
    analyze_stbc_diversity()