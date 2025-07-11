import numpy as np

def matrix_to_str(m):
    """Formats a numpy array into a string for printing."""
    return '\n'.join([' '.join(map(str, row)) for row in m])

def check_matrices():
    """
    Analyzes which of the given matrices are in the set P.
    """
    # Define the matrices
    A = np.array([[0., 0.], [0., 0.]])
    B = np.array([[6., 4.], [3., 7.]])
    C = np.array([[1., -0.5], [-0.5, 1.]])
    D = np.array([[np.pi, 1.], [1., np.pi**2]])
    E = np.array([[1., np.pi], [np.pi, 1.]])
    F = np.array([[42., 0.], [0., 0.]])

    matrices_in_P = []

    print("Analyzing the properties of matrices in P = ConvexHull(v*v^T).")
    print("Any matrix M in P must be: 1. Symmetric, 2. Positive Semidefinite, 3. Have Trace(M) >= 1.")
    print("-" * 40)

    # --- Matrix A ---
    print("--- Checking Matrix A ---\n" + matrix_to_str(A))
    trace_A = np.trace(A)
    print(f"Trace(A) = {trace_A}")
    print("The trace of any matrix in P must be >= 1. Since Trace(A) is 0, A is NOT in P.")
    print("-" * 40)

    # --- Matrix B ---
    print("--- Checking Matrix B ---\n" + matrix_to_str(B))
    is_symmetric_B = np.all(B == B.T)
    print(f"Is B symmetric? {is_symmetric_B}")
    print("A matrix in P must be symmetric. B is not symmetric, so B is NOT in P.")
    print("-" * 40)

    # --- Matrix E ---
    print("--- Checking Matrix E ---\n" + matrix_to_str(E))
    is_psd_E = np.all(np.linalg.eigvalsh(E) >= -1e-9)
    print(f"Is E Positive Semidefinite? {is_psd_E} (det(E) = {np.linalg.det(E):.4f})")
    print("A matrix in P must be Positive Semidefinite. The determinant of E is negative, so it is not PSD. Thus, E is NOT in P.")
    print("-" * 40)

    # --- Matrix C ---
    print("--- Checking Matrix C ---\n" + matrix_to_str(C))
    print("C satisfies the necessary conditions (Symmetric, PSD, Trace >= 1).")
    print("We check if it can be expressed as a convex combination of generators.")
    v1 = np.array([[1.], [1.]])
    S1 = v1 @ v1.T
    v2 = np.array([[1.], [-1.]])
    S2 = v2 @ v2.T
    lambda1 = 1/4
    lambda2 = 3/4
    combination_C = lambda1 * S1 + lambda2 * S2
    print(f"We propose the combination: C = {lambda1} * S(1,1) + {lambda2} * S(1,-1)")
    print(f"Equation: C = {lambda1} * [[{S1[0,0]}, {S1[0,1]}], [{S1[1,0]}, {S1[1,1]}]] + {lambda2} * [[{S2[0,0]}, {S2[0,1]}], [{S2[1,0]}, {S2[1,1]}]]")
    print("Resulting matrix:\n" + matrix_to_str(combination_C))
    if np.allclose(C, combination_C):
        print("The combination is correct. Thus, C is in P.")
        matrices_in_P.append('C')
    print("-" * 40)
    
    # --- Matrix F ---
    print("--- Checking Matrix F ---\n" + matrix_to_str(F))
    print("F satisfies the necessary conditions (Symmetric, PSD, Trace >= 1).")
    print("For a matrix in P, if the (2,2) entry is 0, all generating vectors v=(a,b) must have b=0.")
    print("So F must be a convex combination of S(a,0) = [[a^2, 0], [0, 0]].")
    print("We need to express the (1,1) entry, 42, as a convex combination of non-zero integer squares.")
    print("We can use the squares 36=6^2 and 49=7^2, which bracket 42.")
    v1 = np.array([[6.], [0.]])
    S1 = v1 @ v1.T
    v2 = np.array([[7.], [0.]])
    S2 = v2 @ v2.T
    lambda1 = 7/13
    lambda2 = 6/13
    combination_F = lambda1 * S1 + lambda2 * S2
    print(f"We propose the combination: F = {lambda1:.4f} * S(6,0) + {lambda2:.4f} * S(7,0)")
    print(f"Equation: F = {lambda1:.4f} * [[{S1[0,0]}, {S1[0,1]}], [{S1[1,0]}, {S1[1,1]}]] + {lambda2:.4f} * [[{S2[0,0]}, {S2[0,1]}], [{S2[1,0]}, {S2[1,1]}]]")
    print("Resulting matrix:\n" + matrix_to_str(combination_F))
    if np.allclose(F, combination_F):
        print("The combination is correct. Thus, F is in P.")
        matrices_in_P.append('F')
    print("-" * 40)

    # --- Matrix D ---
    print("--- Checking Matrix D ---\n" + matrix_to_str(D))
    print("D satisfies the necessary conditions (Symmetric, PSD, Trace > 1).")
    print("However, its entries involve the transcendental number pi.")
    print("If D were in P, the vector of its unique entries (pi, 1, pi^2) would be a convex combination of integer vectors (a^2, ab, b^2).")
    print("This implies that pi, 1, and pi^2 would be linearly dependent over the field of rational numbers.")
    print("This would mean pi is a root of a polynomial with rational coefficients (e.g., q_1*pi^2 + q_2*pi + q_3 = 0).")
    print("This contradicts the known mathematical fact that pi is a transcendental number.")
    print("Therefore, D is NOT in P.")
    print("-" * 40)

    # --- Final Answer ---
    print("The list of matrices contained in P is:")
    final_answer = f"[{','.join(sorted(matrices_in_P))}]"
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    check_matrices()