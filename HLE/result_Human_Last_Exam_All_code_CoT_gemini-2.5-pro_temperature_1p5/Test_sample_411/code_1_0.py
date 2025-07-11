import sympy as sp

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three 2x2 space-time block codes.
    """
    # Define symbolic variables for the differences between symbols.
    # d1 = x1 - e1, d2 = x2 - e2. These are complex numbers.
    # We represent d1 as dr1 + I*di1 and d2 as dr2 + I*di2.
    dr1, di1, dr2, di2 = sp.symbols('delta_r1 delta_i1 delta_r2 delta_i2', real=True)
    d1 = dr1 + sp.I * di1
    d2 = dr2 + sp.I * di2

    print("Diversity Analysis of Space-Time Codes")
    print("========================================")
    print("For a 2x1 MISO system, the maximum possible diversity order is 2.")
    print("Diversity is determined by the rank of the codeword difference matrix Delta_S.\n")

    # --- Analysis of Code S_a ---
    print("--- Analysis of Code S_a ---")
    delta_S_a = sp.Matrix([[d1, d2], [d2, d1]])
    det_S_a = sp.det(delta_S_a)
    print("Difference Matrix Delta_S_a:")
    sp.pprint(delta_S_a)
    print(f"\nDeterminant det(Delta_S_a) = {det_S_a}")
    print("The determinant can be zero if d1^2 - d2^2 = 0, which happens if d1 = d2 or d1 = -d2.")
    print("Let's test with a numerical example where d1 = d2.")
    d1_val = 1 + 2j
    d2_val = 1 + 2j
    det_val_a = d1_val**2 - d2_val**2
    print(f"Let d1 = {d1_val} and d2 = {d2_val}. Then, det = ({d1_val})^2 - ({d2_val})^2 = {d1_val**2} - {d2_val**2} = {det_val_a}")
    print("Since the determinant is zero for distinct codewords, the matrix is not always full rank.")
    print("Diversity order for S_a is 1.\n")

    # --- Analysis of Code S_b ---
    print("--- Analysis of Code S_b ---")
    delta_S_b = sp.Matrix([[d1, d2], [d2, sp.conjugate(d1)]])
    det_S_b = sp.simplify(sp.det(delta_S_b))
    print("Difference Matrix Delta_S_b:")
    sp.pprint(delta_S_b)
    print(f"\nDeterminant det(Delta_S_b) = {det_S_b}")
    print("The determinant can be zero if |d1|^2 - d2^2 = 0.")
    print("This can happen if d2 is real and d2^2 = |d1|^2. Such values can be constructed from M-QAM symbols.")
    print("Let's test with a numerical example: d1 = 2j and d2 = 2.")
    d1_val_b = 2j
    d2_val_b = 2
    det_val_b = abs(d1_val_b)**2 - d2_val_b**2
    print(f"Let d1 = {d1_val_b} and d2 = {d2_val_b}. Then, det = |{d1_val_b}|^2 - ({d2_val_b})^2 = {abs(d1_val_b)**2} - {d2_val_b**2} = {det_val_b}")
    print("Since the determinant is zero for distinct codewords, the matrix is not always full rank.")
    print("Diversity order for S_b is 1.\n")
    
    # --- Analysis of Code S_c ---
    print("--- Analysis of Code S_c ---")
    delta_S_c = sp.Matrix([[-sp.conjugate(d1), d2], [-sp.conjugate(d2), -d1]])
    det_S_c = sp.simplify(sp.det(delta_S_c))
    print("Difference Matrix Delta_S_c:")
    sp.pprint(delta_S_c)
    print(f"\nDeterminant det(Delta_S_c) = {det_S_c}")
    print("The determinant |d1|^2 + |d2|^2 is a sum of non-negative terms.")
    print("It is zero only if both d1 = 0 and d2 = 0, which means the codewords are identical, not distinct.")
    print("For any pair of distinct codewords, the determinant is always greater than zero.")
    print("Thus, the difference matrix is always full rank.")
    print("Diversity order for S_c is 2.\n")

    # --- Conclusion ---
    print("--- Final Conclusion ---")
    print("(a) The diversity order for each code is:")
    print("    - S_a: 1")
    print("    - S_b: 1")
    print("    - S_c: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("Assuming 'directivity' refers to the diversity order, code S_c provides the maximum directivity (2).")


if __name__ == '__main__':
    analyze_stbc_diversity()
