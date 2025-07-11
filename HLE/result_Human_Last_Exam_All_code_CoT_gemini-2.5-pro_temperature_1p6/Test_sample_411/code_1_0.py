import sympy

def analyze_space_time_codes():
    """
    Analyzes the diversity order of three space-time codes for a 2x1 MISO system
    by checking the rank criterion on the codeword difference matrix.
    """
    # --- Introduction to the method ---
    print("--- Diversity Order Analysis ---")
    print("The diversity order of a space-time code is determined by the rank of the codeword difference matrix,")
    print("ΔS = S_i - S_j. For a 2x1 MISO system, full diversity (order 2) is achieved if ΔS is always")
    print("full rank (i.e., det(ΔS) != 0) for any two distinct codewords S_i and S_j.\n")

    # Define symbolic variables for the symbol differences dx1 = x1 - x'_1 and dx2 = x2 - x'_2
    dx1 = sympy.Symbol('Δx1', complex=True)
    dx2 = sympy.Symbol('Δx2', complex=True)

    # --- (a) Analysis of Code Sa ---
    print("--- Analysis of Code S_a ---")
    print("S_a = [[x1, x2], [x2, x1]]")
    Delta_Sa = sympy.Matrix([[dx1, dx2], [dx2, dx1]])
    print("The difference matrix is ΔS_a = [[Δx1, Δx2], [Δx2, Δx1]]")

    det_Sa = Delta_Sa.det()
    print(f"\nThe determinant is det(ΔS_a) = {det_Sa}")
    print("This determinant becomes zero if Δx1² = Δx2². This can occur for non-zero differences.")
    print("For example, if Δx1 = 2 and Δx2 = 2 (a possible difference for BPSK symbols):")
    val_dx1_a = 2
    val_dx2_a = 2
    det_val_a = val_dx1_a**2 - val_dx2_a**2
    print(f"det(ΔS_a) = {val_dx1_a}² - {val_dx2_a}² = {val_dx1_a**2} - {val_dx2_a**2} = {det_val_a}")
    print("Since the determinant can be zero for distinct codewords, the diversity order is 1.\n")

    # --- (a) Analysis of Code Sb ---
    print("--- Analysis of Code S_b ---")
    print("S_b = [[x1, x2], [x2, x1*]]")
    print("The difference matrix is ΔS_b = [[Δx1, Δx2], [Δx2, Δx1*]]")

    # Note: sympy.conjugate() represents the complex conjugate.
    det_Sb = sympy.det(sympy.Matrix([[dx1, dx2], [dx2, dx1.conjugate()]]))
    print(f"\nThe determinant is det(ΔS_b) = {det_Sb}, which simplifies to |Δx1|² - Δx2².")
    print("This determinant can also be zero for non-zero differences.")
    print("For example, consider differences from a QPSK system. Let Δx1 = 2j and Δx2 = 2:")
    val_dx1_b = 2j
    val_dx2_b = 2
    v1_re, v1_im = int(val_dx1_b.real), int(val_dx1_b.imag)
    abs1_sq = abs(val_dx1_b)**2
    det_val_b = abs1_sq - val_dx2_b**2
    print(f"det(ΔS_b) = |{val_dx1_b}|² - ({val_dx2_b})² = ({v1_re}² + {v1_im}²) - {val_dx2_b}² = {int(abs1_sq)} - {val_dx2_b**2} = {int(det_val_b)}")
    print("Since the determinant can be zero for distinct codewords, the diversity order is 1.\n")

    # --- (a) Analysis of Code Sc ---
    print("--- Analysis of Code S_c ---")
    print("S_c = [[-x1*, x2], [-x2*, -x1]]")
    print("The difference matrix is ΔS_c = [[-Δx1*, Δx2], [-Δx2*, -Δx1]]")
    
    det_Sc = sympy.simplify(sympy.det(sympy.Matrix([[-dx1.conjugate(), dx2], [-dx2.conjugate(), -dx1]])))
    print(f"\nThe determinant is det(ΔS_c) = {det_Sc}, which simplifies to |Δx1|² + |Δx2|².")
    print("This expression is the sum of two non-negative terms. It is zero if and only if")
    print("both terms are zero, which means Δx1 = 0 and Δx2 = 0.")
    print("This only happens if the codewords are identical, not distinct.")
    print("For any distinct codewords, the determinant is always non-zero. For example, if Δx1 = 1+1j, Δx2 = 2-3j:")
    val_dx1_c = 1+1j
    val_dx2_c = 2-3j
    v1_re, v1_im = int(val_dx1_c.real), int(val_dx1_c.imag)
    v2_re, v2_im = int(val_dx2_c.real), int(val_dx2_c.imag)
    abs1_sq = abs(val_dx1_c)**2
    abs2_sq = abs(val_dx2_c)**2
    print(f"det(ΔS_c) = |{val_dx1_c}|² + |{val_dx2_c}|² = ({v1_re}² + {v1_im}²) + ({v2_re}² + ({v2_im})²) = {int(abs1_sq)} + {int(abs2_sq)} = {int(abs1_sq + abs2_sq)}")
    print("Since the determinant is always non-zero, the code achieves full rank. The diversity order is 2.\n")
    
    # --- Final Summary ---
    print("--- Summary and Conclusion ---")
    print("(a) Diversity Orders:")
    print(" - Code S_a: Diversity Order = 1")
    print(" - Code S_b: Diversity Order = 1")
    print(" - Code S_c: Diversity Order = 2 (Full Diversity)")
    print("\n(b) Maximum Directivity:")
    print("Code S_c provides the maximum directivity (diversity).")

if __name__ == '__main__':
    analyze_space_time_codes()