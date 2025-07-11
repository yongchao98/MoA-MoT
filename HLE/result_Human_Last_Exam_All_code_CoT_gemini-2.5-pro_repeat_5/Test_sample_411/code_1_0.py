import sympy

def analyze_stbc_diversity():
    """
    This function analyzes the diversity order of three 2x2 space-time block codes
    and determines which one provides the maximum diversity for a 2x1 MISO system.
    """
    # Introduction to the method
    print("--- Diversity Order Analysis of Space-Time Block Codes ---")
    print("The diversity order of a space-time code is the minimum rank of the difference")
    print("matrix Delta = S - E for any two distinct codewords S and E.")
    print("A code provides full diversity if this rank is equal to the number of transmit")
    print("antennas (N=2 in this case). This is guaranteed if the determinant of Delta")
    print("is non-zero for any non-zero symbol difference vector (delta_1, delta_2).\n")

    # Define symbolic variables for the differences in transmitted symbols
    # d1 = x1 - e1, d2 = x2 - e2
    d1, d2 = sympy.symbols('delta_1, delta_2', complex=True)

    # --- Analysis of Code Sa ---
    print("--- (1) Analysis of Code Sa ---")
    # Delta_a = [ [d1, d2], [d2, d1] ]
    Delta_a = sympy.Matrix([[d1, d2], [d2, d1]])
    det_a = Delta_a.det()
    print("Difference Matrix Delta_a:")
    sympy.pprint(Delta_a)
    print(f"\nDeterminant of Delta_a: det(Delta_a) = {sympy.pretty(det_a)}")
    print("The determinant is zero if delta_1^2 = delta_2^2, which can occur for non-zero deltas")
    print("(e.g., if delta_1 = delta_2). This means the matrix can be rank-deficient.")
    print("Result: The diversity order for Sa is 1.\n")
    
    print("Example of a rank-deficient matrix for Sa:")
    d1_val, d2_val = 1, 1
    example_a = Delta_a.subs({d1: d1_val, d2: d2_val})
    print(f"Let delta_1 = {d1_val} and delta_2 = {d2_val}. The matrix is:")
    sympy.pprint(example_a)
    print("The final equation for the determinant is:")
    print(f"det = ({d1_val}) * ({d1_val}) - ({d2_val}) * ({d2_val}) = {d1_val*d1_val - d2_val*d2_val}")
    print(f"Since the determinant is 0, the rank is {example_a.rank()}.\n")

    # --- Analysis of Code Sb ---
    print("--- (2) Analysis of Code Sb ---")
    # Delta_b = [ [d1, d2], [d2, d1*] ]
    Delta_b = sympy.Matrix([[d1, d2], [d2, sympy.conjugate(d1)]])
    det_b = Delta_b.det()
    print("Difference Matrix Delta_b:")
    sympy.pprint(Delta_b)
    print(f"\nDeterminant of Delta_b: det(Delta_b) = {sympy.pretty(det_b)}")
    print("The determinant is zero if |delta_1|^2 = delta_2^2. This can occur for non-zero deltas")
    print("(e.g., if delta_1 is complex and delta_2 is real with equal magnitude).")
    print("Result: The diversity order for Sb is 1.\n")

    print("Example of a rank-deficient matrix for Sb:")
    d1_val, d2_val = 3 + 4j, 5
    example_b = Delta_b.subs({d1: d1_val, d2: d2_val})
    print(f"Let delta_1 = {d1_val} and delta_2 = {d2_val}. The matrix is:")
    sympy.pprint(example_b)
    det_b_val = d1_val * sympy.conjugate(d1_val) - d2_val**2
    print("The final equation for the determinant is:")
    print(f"det = ({d1_val}) * ({sympy.conjugate(d1_val)}) - ({d2_val}) * ({d2_val}) = {sympy.N(det_b_val)}")
    print(f"Since the determinant is 0, the rank is {example_b.rank()}.\n")

    # --- Analysis of Code Sc ---
    print("--- (3) Analysis of Code Sc ---")
    # Delta_c = [ [-d1*, d2], [-d2*, -d1] ]
    Delta_c = sympy.Matrix([[-sympy.conjugate(d1), d2], [-sympy.conjugate(d2), -d1]])
    det_c = sympy.expand(Delta_c.det())
    print("Difference Matrix Delta_c:")
    sympy.pprint(Delta_c)
    print(f"\nDeterminant of Delta_c: det(Delta_c) = {sympy.pretty(det_c)} = |delta_1|^2 + |delta_2|^2")
    print("This determinant is a sum of squares of magnitudes. It is zero if and only if")
    print("both delta_1 = 0 and delta_2 = 0. For any distinct codewords (non-zero deltas),")
    print("the determinant is non-zero, and the matrix is full rank.")
    print("Result: The diversity order for Sc is 2 (full diversity).\n")
    print("The final equation for the determinant of Sc shows it is always non-negative:")
    print("det = (-conj(delta_1))*(-delta_1) - (delta_2)*(-conj(delta_2)) = |delta_1|^2 + |delta_2|^2\n")

    # --- Final Summary ---
    print("--- Summary of Results ---")
    print("(a) What is the diversity order for each code?")
    print("    - Diversity Order of Sa: 1")
    print("    - Diversity Order of Sb: 1")
    print("    - Diversity Order of Sc: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("    The maximum possible diversity order for this system is 2.")
    print("    Code Sc is the only code that achieves this maximum diversity.")
    print("    Therefore, Sc provides the maximum 'directivity' (diversity gain).")

if __name__ == '__main__':
    analyze_stbc_diversity()

<<<S_c>>>