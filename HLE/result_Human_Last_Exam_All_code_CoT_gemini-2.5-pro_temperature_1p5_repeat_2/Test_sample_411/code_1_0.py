def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three 2x2 space-time block codes.
    """
    print("--- Analysis of Space-Time Block Code Diversity Order ---")
    print("\nThe diversity order is determined by the minimum rank of the codeword difference matrix Delta_S.")
    print("For a 2x2 code, a diversity order of 2 (full diversity) is achieved if det(Delta_S) is non-zero for any non-zero error.\n")
    print("Let d1 = x1 - e1 and d2 = x2 - e2. The error is non-zero if (d1, d2) != (0, 0).\n")

    # --- Code Sa ---
    print("--- (1) Analysis of Code Sa ---")
    print("S_a = [[x1, x2],")
    print("       [x2, x1]]")
    print("The difference matrix is Delta_S_a = [[d1, d2], [d2, d1]].")
    print("The determinant is det(Delta_S_a) = (d1 * d1) - (d2 * d2) = d1^2 - d2^2.")
    print("This determinant can be zero for non-zero errors. For example, if d1 = d2 = 1, then d1^2 - d2^2 = 1 - 1 = 0.")
    print("When the determinant is zero, the matrix is singular and its rank is 1.")
    print("Result: The minimum rank is 1, so the diversity order for S_a is 1.\n")

    # --- Code Sb ---
    print("--- (2) Analysis of Code Sb ---")
    print("S_b = [[x1,   x2],")
    print("       [x2, x1*]]  (where x1* is the complex conjugate of x1)")
    print("The difference matrix is Delta_S_b = [[d1, d2], [d2, d1*]].")
    print("The determinant is det(Delta_S_b) = (d1 * d1*) - (d2 * d2) = |d1|^2 - d2^2.")
    print("This determinant can also be zero for non-zero errors. For example, if d1 = 2 and d2 = 2 (symbols from a QAM constellation can result in real-valued differences), then |d1|^2 - d2^2 = |2|^2 - 2^2 = 4 - 4 = 0.")
    print("When the determinant is zero, the matrix rank is 1.")
    print("Result: The minimum rank is 1, so the diversity order for S_b is 1.\n")

    # --- Code Sc ---
    print("--- (3) Analysis of Code Sc ---")
    print("S_c = [[-x1*, x2],")
    print("       [-x2*, -x1]]")
    print("The difference matrix is Delta_S_c = [[-d1*,  d2], [-d2*, -d1]].")
    print("The determinant is det(Delta_S_c) = (-d1*) * (-d1) - (d2) * (-d2*) = |d1|^2 + |d2|^2.")
    print("This determinant is a sum of two non-negative terms. It can only be zero if |d1|^2 = 0 AND |d2|^2 = 0, which means d1 = 0 AND d2 = 0.")
    print("This corresponds to a zero error, so for any non-zero error, the determinant is strictly positive.")
    print("Therefore, the matrix is always full rank (rank 2) for any non-zero error event.")
    print("Result: The minimum rank is 2, so the diversity order for S_c is 2 (full diversity).\n")

    # --- Summary and Conclusion ---
    print("--- Summary and Conclusion ---")
    print("(a) What is the diversity order for each code?")
    print("    - Diversity order of S_a: 1")
    print("    - Diversity order of S_b: 1")
    print("    - Diversity order of S_c: 2\n")

    print("(b) Which code provides the maximum directivity?")
    print("Assuming 'directivity' refers to diversity gain, the code with the highest diversity order provides the best performance against fading.")
    print("Code S_c achieves the maximum possible diversity order of 2.")
    print("Conclusion: Code S_c provides the maximum directivity.\n")

if __name__ == '__main__':
    analyze_stbc_diversity()
