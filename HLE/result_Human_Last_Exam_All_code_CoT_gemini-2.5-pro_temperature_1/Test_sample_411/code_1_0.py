import sympy as sp

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three space-time block codes (STBCs)
    and identifies the one with the maximum directivity.
    """
    print("This script analyzes the diversity order of three STBCs for a 2x1 MISO system.")
    print("The diversity order is the minimum rank of the codeword difference matrix, Delta_S.")
    print("We analyze this by checking if the determinant of Delta_S can be zero for")
    print("non-zero symbol differences (dx1, dx2), where dx1=x1-e1 and dx2=x2-e2.")
    print("-" * 70)

    # Define symbolic variables for the differences of the symbols.
    # We use real and imaginary parts for clear handling of complex conjugates.
    dx1_re, dx1_im = sp.symbols('dx1_re dx1_im', real=True)
    dx2_re, dx2_im = sp.symbols('dx2_re dx2_im', real=True)
    dx1 = dx1_re + sp.I * dx1_im
    dx2 = dx2_re + sp.I * dx2_im
    dx1_c = sp.conjugate(dx1)
    dx2_c = sp.conjugate(dx2)

    # --- Part (a): Diversity Order Analysis ---
    print("\n(a) Diversity Order Analysis")

    # --- Analysis for Code Sa ---
    print("\n--- Code S_a: S = [[x1, x2], [x2, x1]] ---")
    # Difference matrix: [[dx1, dx2], [dx2, dx1]]
    det_Sa = dx1**2 - dx2**2
    print("The determinant of the difference matrix Delta_S_a is: det = dx1**2 - dx2**2")
    print("In terms of real and imaginary parts, the equation for the determinant being zero is:")
    print(f"{sp.re(det_Sa).expand()} = 0  (real part)")
    print(f"{sp.im(det_Sa).expand()} = 0  (imaginary part)")
    print("\nThe determinant is zero if dx1 = dx2 or dx1 = -dx2. This is possible for")
    print("distinct symbols from a QAM constellation (e.g., dx1=2, dx2=2).")
    print("Since the determinant can be zero, the minimum rank is 1.")
    diversity_a = 1
    print(f"The diversity order for code S_a is: {diversity_a}")
    print("-" * 70)

    # --- Analysis for Code Sb ---
    print("\n--- Code S_b: S = [[x1, x2], [x2, x1*]] ---")
    # Difference matrix: [[dx1, dx2], [dx2, dx1*]]
    det_Sb = dx1_c * dx1 - dx2**2
    print("The determinant of the difference matrix Delta_S_b is: det = |dx1|^2 - dx2**2")
    print("In terms of real and imaginary parts, the equation for the determinant being zero is:")
    print(f"{sp.re(det_Sb).expand()} = 0  (real part)")
    print(f"{sp.im(det_Sb).expand()} = 0  (imaginary part)")
    print("\nThe determinant is zero if |dx1|^2 = dx2**2. The imaginary part of this equation requires")
    print("dx2 to be purely real. The real part requires |dx1|^2 = (real_dx2)^2.")
    print("This is possible (e.g., dx1=2, dx2=2).")
    print("Since the determinant can be zero, the minimum rank is 1.")
    diversity_b = 1
    print(f"The diversity order for code S_b is: {diversity_b}")
    print("-" * 70)

    # --- Analysis for Code Sc ---
    print("\n--- Code S_c: S = [[-x1*, x2], [-x2*, -x1]] ---")
    # Difference matrix: [[-dx1*, dx2], [-dx2*, -dx1]]
    det_Sc = (-dx1_c) * (-dx1) - (dx2) * (-dx2_c)
    print("The determinant of the difference matrix Delta_S_c is: det = |dx1|^2 + |dx2|^2")
    print("In terms of real and imaginary parts, the equation for the determinant is:")
    print(f"det = {det_Sc.expand()}")
    print("\nThe determinant is a sum of non-negative terms. It is zero only if |dx1|^2=0 and |dx2|^2=0,")
    print("which means dx1=0 and dx2=0. This only happens for identical codewords, which are excluded.")
    print("Therefore, the determinant is always non-zero for distinct codewords, so the matrix always has full rank (2).")
    diversity_c = 2
    print(f"The diversity order for code S_c is: {diversity_c}")
    print("-" * 70)

    # --- Part (b): Maximum Directivity ---
    print("\n(b) Maximum Directivity")
    print("Maximum directivity, in the context of combating fading, corresponds to the maximum diversity order.")
    print(f"The maximum possible diversity order for a 2x1 system is NL = 2*1 = 2.")
    print(f"Code S_c is the only code that achieves this maximum diversity order of {diversity_c}.")
    print("Therefore, code S_c provides the maximum directivity.")

if __name__ == '__main__':
    analyze_stbc_diversity()
