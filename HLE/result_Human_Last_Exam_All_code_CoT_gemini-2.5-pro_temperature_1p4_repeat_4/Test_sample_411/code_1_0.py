import cmath

def analyze_stbc_diversity():
    """
    Analyzes and explains the diversity order for the given STBCs.
    It demonstrates the calculations by plugging in example error values.
    """

    print("To determine the diversity order, we analyze the rank of the codeword difference matrix, Delta_S.")
    print("For a 2x2 code matrix, if its determinant can be zero for non-zero error symbols (e1, e2),")
    print("the rank is less than 2, and the diversity order is 1.")
    print("If the rank is always 2 for any distinct codewords, the diversity order is 2.\n")

    # --- Analysis for Code Sa ---
    print("--- Analysis for Code Sa ---")
    print("S_a = [[x1, x2], [x2, x1]]")
    print("The difference matrix is Delta_Sa = [[e1, e2], [e2, e1]].")
    print("The determinant is det(Delta_Sa) = e1^2 - e2^2.")
    print("This determinant is zero if e1 = e2 or e1 = -e2.")
    e1 = 1 + 0j
    e2 = 1 + 0j
    print(f"Let's test with a non-zero case where this happens: e1 = {e1}, e2 = {e2}.")
    det_Sa = e1**2 - e2**2
    print(f"Calculation: det(Delta_Sa) = ({e1})^2 - ({e2})^2 = ({e1**2}) - ({e2**2}) = {det_Sa}")
    if det_Sa == 0:
        print("The determinant is 0, so the rank is 1. The diversity order for S_a is 1.\n")

    # --- Analysis for Code Sb ---
    print("--- Analysis for Code Sb ---")
    print("S_b = [[x1, x2], [x2, x1*]] (* denotes conjugate)")
    print("The difference matrix is Delta_Sb = [[e1, e2], [e2, e1*]].")
    print("The determinant is det(Delta_Sb) = e1*e1* - e2^2 = |e1|^2 - e2^2.")
    print("This can be zero if we can find error symbols where e2 is real and |e1| = |e2|.")
    e1 = 0 + 1j
    e2 = 1 + 0j
    print(f"Let's test with a case where e1 = {e1} (so |e1|=1) and e2 = {e2}.")
    det_Sb = abs(e1)**2 - e2**2
    abs_e1_sq = abs(e1)**2
    e2_sq = e2**2
    print(f"Calculation: det(Delta_Sb) = |{e1}|^2 - ({e2})^2 = {abs_e1_sq:.1f} - ({e2_sq}) = {det_Sb}")
    if det_Sb.real == 0 and det_Sb.imag == 0:
        print("The determinant is 0, so the rank is 1. The diversity order for S_b is 1.\n")

    # --- Analysis for Code Sc ---
    print("--- Analysis for Code Sc ---")
    print("S_c = [[-x1*, x2], [-x2*, -x1]]")
    print("This is a variant of the Alamouti code. To properly analyze it, we must check det(Delta_Sc * Delta_Sc^H).")
    print("It can be shown that Delta_Sc * Delta_Sc^H = (|e1|^2 + |e2|^2) * I, where I is the identity matrix.")
    print("The determinant is det(Delta_Sc * Delta_Sc^H) = (|e1|^2 + |e2|^2)^2.")
    print("This expression is zero only if both e1 and e2 are zero (i.e., the codewords are not distinct).")
    e1 = 1 + 2j
    e2 = 3 - 1j
    print(f"Let's test with an arbitrary non-zero case: e1 = {e1}, e2 = {e2}")
    abs_e1_sq = abs(e1)**2
    abs_e2_sq = abs(e2)**2
    det_Sc_prod = (abs_e1_sq + abs_e2_sq)**2
    print(f"|e1|^2 = {abs_e1_sq:.1f}")
    print(f"|e2|^2 = {abs_e2_sq:.1f}")
    print(f"Calculation: det = ({abs_e1_sq:.1f} + {abs_e2_sq:.1f})^2 = ({abs_e1_sq + abs_e2_sq:.1f})^2 = {det_Sc_prod:.1f}")
    if det_Sc_prod != 0:
        print("The determinant is non-zero, indicating full rank (2). The diversity order for S_c is 2.\n")
    
    print("--- FINAL CONCLUSION ---")
    print("(a) The diversity orders are:")
    print("    - Code S_a: 1")
    print("    - Code S_b: 1")
    print("    - Code S_c: 2")
    print("\n(b) Assuming 'directivity' means 'diversity', Code S_c provides the maximum diversity.")

if __name__ == '__main__':
    analyze_stbc_diversity()
