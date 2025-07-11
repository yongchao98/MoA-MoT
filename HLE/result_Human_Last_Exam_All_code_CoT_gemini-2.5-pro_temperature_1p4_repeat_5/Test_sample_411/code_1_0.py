import sympy as sp

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three space-time codes Sa, Sb, and Sc.
    """
    # Step 1: Define symbolic variables for symbol differences.
    # dx1 = x1 - x1', dx2 = x2 - x2'. These are complex numbers.
    # A non-zero difference vector [dx1, dx2] corresponds to two distinct codewords.
    dx1, dx2 = sp.symbols('dx1 dx2', complex=True)

    # Step 2: Define the codeword difference matrices (dS).
    # The diversity order is the minimum rank of dS over all non-zero [dx1, dx2].
    
    # Code Sa
    dS_a = sp.Matrix([
        [dx1, dx2],
        [dx2, dx1]
    ])

    # Code Sb
    dS_b = sp.Matrix([
        [dx1, dx2],
        [dx2, sp.conjugate(dx1)]
    ])

    # Code Sc
    dS_c = sp.Matrix([
        [-sp.conjugate(dx1), dx2],
        [-sp.conjugate(dx2), -dx1]
    ])

    # Step 3: Calculate the determinant of each difference matrix.
    # If det(dS) can be zero for non-zero dx1 or dx2, the code is not full diversity.
    det_a = dS_a.det()
    det_b = dS_b.det()
    det_c = dS_c.det()
    
    print("--- Analysis of Space-Time Code Diversity ---")
    print("We check if the determinant of the codeword difference matrix dS can be zero for non-zero inputs.")
    print("If det(dS) is always non-zero, the code achieves full diversity (order 2).\n")

    # --- Analysis for Code Sa ---
    print("--- Code Sa ---")
    print("Matrix dS_a:")
    print(sp.pretty(dS_a))
    print(f"\ndet(dS_a) = {det_a}")
    print("This determinant equals zero if dx1^2 = dx2^2, which can happen for non-zero inputs.")
    # Example where det(dS_a) is zero
    ex_dx1_a = 2
    ex_dx2_a = 2
    det_a_val = det_a.subs({dx1: ex_dx1_a, dx2: ex_dx2_a})
    print(f"Example: Let dx1 = {ex_dx1_a} and dx2 = {ex_dx2_a}. Both are non-zero.")
    print(f"det(dS_a) = ({ex_dx1_a})^2 - ({ex_dx2_a})^2 = {ex_dx1_a**2} - {ex_dx2_a**2} = {det_a_val}")
    print("Since the determinant can be zero, the rank of dS_a can be 1.")
    print("Diversity order of Sa is 1.\n")

    # --- Analysis for Code Sb ---
    print("--- Code Sb ---")
    print("Matrix dS_b:")
    print(sp.pretty(dS_b))
    print(f"\ndet(dS_b) = {det_b}")
    print("This determinant equals zero if |dx1|^2 = dx2^2. This can also happen for non-zero inputs.")
    # Example where det(dS_b) is zero
    ex_dx1_b = 2.0
    ex_dx2_b = 2.0
    det_b_val = (abs(ex_dx1_b)**2 - ex_dx2_b**2)
    print(f"Example: Let dx1 = {ex_dx1_b} and dx2 = {ex_dx2_b}. Both are non-zero.")
    print(f"det(dS_b) = |{ex_dx1_b}|^2 - ({ex_dx2_b})^2 = {abs(ex_dx1_b)**2} - {ex_dx2_b**2} = {det_b_val}")
    print("Since the determinant can be zero, the rank of dS_b can be 1.")
    print("Diversity order of Sb is 1.\n")

    # --- Analysis for Code Sc ---
    print("--- Code Sc ---")
    print("Matrix dS_c:")
    print(sp.pretty(dS_c))
    det_c_expanded = sp.expand(det_c)
    print(f"\ndet(dS_c) = {det_c_expanded}")
    print("This expression simplifies to |dx1|^2 + |dx2|^2.")
    print("The sum of two non-negative terms |dx1|^2 and |dx2|^2 is zero only if dx1 = 0 and dx2 = 0.")
    print("For any distinct codewords (non-zero dx1 or dx2), the determinant is always positive.")
    print("Therefore, the matrix dS_c is always full rank (rank 2).")
    print("Diversity order of Sc is 2 (Full Diversity).\n")

    # Step 4: Final conclusions for the problem.
    print("--- Final Answers ---")
    print("(a) What is the diversity order for each code?")
    print("    - Diversity order of Sa: 1")
    print("    - Diversity order of Sb: 1")
    print("    - Diversity order of Sc: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("    Code Sc provides the maximum directivity because it achieves the highest possible diversity order (full diversity).")

analyze_stbc_diversity()