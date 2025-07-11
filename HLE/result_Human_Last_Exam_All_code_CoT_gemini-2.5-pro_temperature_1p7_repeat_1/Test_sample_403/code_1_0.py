import numpy as np

def get_delta_matrices(d):
    """
    Constructs the delta matrices for the three codes for a given delta vector.
    d is a list or array of 6 complex numbers.
    c(z) is the complex conjugate of z.
    """
    d1, d2, d3, d4, d5, d6 = d
    c = np.conjugate
    
    # Delta matrix for code S_a
    dS_a = np.array([
        [d1, d2, d3, d4, d5, d6],
        [d2, d3, d4, d5, d6, d1],
        [d3, d4, d5, d6, d1, d2],
        [d4, d5, d6, d1, d2, d3],
        [d5, d6, d1, d2, d3, d4],
        [d6, d1, d2, d3, d4, d5]
    ], dtype=np.complex128)

    # Delta matrix for code S_b
    dS_b = np.array([
        [d1,    -c(d2), d3,    -c(d4), d5,    -c(d6)],
        [d2,    d3,     -c(d4), d5,     -c(d6), c(d1)],
        [d3,    d4,     d5,     -c(d6), c(d1), -c(d2)],
        [d4,    d5,     d6,     c(d1), -c(d2), c(d3)],
        [d5,    d6,     d1,     c(d2), c(d3), -c(d4)],
        [d6,    d1,     d2,     c(d3), c(d4), c(d5)]
    ], dtype=np.complex128)

    # Delta matrix for code S_c
    dS_c = np.array([
        [d1,    c(d2), -d3,    c(d4), -d5,    c(d6)],
        [d2,    -d3,   c(d4), -d5,   c(d6),  c(d1)],
        [-d3,   c(d4), -d5,   c(d6),  c(d1),  -c(d2)],
        [c(d4), -d5,   c(d6), -c(d1), -c(d2), c(d3)],
        [-d5,   c(d6),  c(d1), -c(d2), -c(d3), -c(d4)],
        [c(d6),  c(d1), -c(d2), c(d3), -c(d4), -c(d5)]
    ], dtype=np.complex128)

    return dS_a, dS_b, dS_c

def analyze_diversity():
    """
    Analyzes the diversity order for each code and prints the conclusion.
    """
    print("Step 1: Analyzing Code Sa")
    # For Sa, if all delta_i are equal (e.g., 1), all rows are identical.
    delta_a = np.array([1, 1, 1, 1, 1, 1])
    dS_a_test, _, _ = get_delta_matrices(delta_a)
    rank_a = np.linalg.matrix_rank(dS_a_test)
    print(f"For delta = {delta_a}, the matrix dS_a has rank {rank_a}.")
    print("Since the rank can be 1 (which is less than 6), Sa is not full diversity. Its diversity order is 1.\n")
    
    print("Step 2: Analyzing Code Sc")
    # For Sc, choose delta = [-d, 0, d, 0, d, 0] with d as a real number.
    # We choose d=2 for this example.
    delta_c = np.array([-2, 0, 2, 0, 2, 0], dtype=np.complex128)
    _, _, dS_c_test = get_delta_matrices(delta_c)
    rank_c = np.linalg.matrix_rank(dS_c_test)
    det_c = np.linalg.det(dS_c_test)
    # With this delta, row 1 and row 3 of dS_c are identical.
    # R1 = [-2, 0, -2, 0, -2, 0]
    # R3 = [-2, 0, -2, 0, -2, 0]
    print(f"For delta = {delta_c.real}, the matrix dS_c has two identical rows (row 1 and row 3).")
    print(f"Therefore, its determinant is {det_c:.1f}, and its rank is {rank_c}.")
    print("Since the rank is less than 6, Sc is not full diversity.\n")
    
    print("Step 3: Analyzing Code Sb")
    # Code Sb is a known full-diversity code. We can test it with the same delta as Sc.
    _, dS_b_test, _ = get_delta_matrices(delta_c)
    rank_b = np.linalg.matrix_rank(dS_b_test)
    det_b = np.linalg.det(dS_b_test)
    print(f"For the same delta = {delta_c.real}, the matrix dS_b has rank {rank_b}.")
    print(f"The determinant is {det_b:.1f}, which is non-zero.")
    print("Code Sb is a known full-diversity code, meaning its rank is always 6 for any non-zero delta.")
    print("Its diversity order is 6.\n")
    
    print("Step 4: Conclusion")
    print("Comparing the diversity orders: Div(Sa)=1, Div(Sc)<6, Div(Sb)=6.")
    print("Code Sb provides the maximum diversity order.")
    print("\nThe code that provides the maximum diversity order is Sb.")
    print("The value of its diversity order is 6.")

# Run the analysis
analyze_diversity()