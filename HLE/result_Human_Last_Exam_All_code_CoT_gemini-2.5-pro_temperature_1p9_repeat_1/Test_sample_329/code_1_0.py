import numpy as np

def is_symmetric(M):
    return np.all(M == M.T)

def is_psd(M):
    if not is_symmetric(M):
        return False
    # For a 2x2 symmetric matrix, PSD check is simpler.
    # But using eigenvalues is general.
    # Note: adding a small tolerance for numerical stability, though not strictly needed here
    # since matrices have integer or exact fraction entries except for D and E.
    eigenvalues = np.linalg.eigvalsh(M)
    return np.all(eigenvalues >= -1e-9)

def check_matrices():
    A = np.array([[0, 0], [0, 0]])
    B = np.array([[6, 4], [3, 7]])
    C = np.array([[1, -1/2], [-1/2, 1]])
    # For D and E, use symbolic or high-precision math conceptually,
    # but for code demonstration, float approximation is fine for some checks.
    pi = np.pi
    D = np.array([[pi, 1], [1, pi**2]])
    E = np.array([[1, pi], [pi, 1]])
    F = np.array([[42, 0], [0, 0]])

    matrices = {'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F}
    in_P = []

    print("--- Analyzing each matrix ---")

    # Matrix A
    print("\nMatrix A:\n", A)
    trace_A = np.trace(A)
    print(f"Trace is {trace_A}. All matrices in P must have trace >= 1.")
    print("Conclusion: A is not in P.\n")

    # Matrix B
    print("\nMatrix B:\n", B)
    if not is_symmetric(B):
        print("Matrix B is not symmetric.")
    print("Conclusion: B is not in P.\n")

    # Matrix C
    print("\nMatrix C:\n", C)
    print(f"Symmetric: {is_symmetric(C)}")
    print(f"PSD: {is_psd(C)}")
    print(f"Trace: {np.trace(C)} >= 1")
    print("C satisfies the necessary conditions. Let's try to construct it.")
    # We try to write C as a convex combination of S_v matrices.
    # Let's try with v1=(1,1) and v2=(1,-1).
    # S1 = [[1, 1], [1, 1]], S2 = [[1, -1], [-1, 1]]
    # C = l1*S1 + l2*S2, with l1+l2=1, l1>=0, l2>=0
    # From off-diagonal: -1/2 = l1*1 + l2*(-1) = l1 - l2
    # Solving l1+l2=1 and l1-l2=-1/2 gives l1=1/4, l2=3/4.
    l1, l2 = 1/4, 3/4
    S1 = np.array([[1, 1], [1, 1]]) # from v=(1,1)
    S2 = np.array([[1, -1], [-1, 1]]) # from v=(1,-1)
    constructed_C = l1 * S1 + l2 * S2
    if np.allclose(C, constructed_C):
        print("C can be constructed as a convex combination:")
        print(f"C = {l1} * {S1.tolist()} + {l2} * {S2.tolist()}")
        print("Conclusion: C is in P.\n")
        in_P.append('C')
    else:
        print("Conclusion: C is not in P.\n")

    # Matrix D
    print("\nMatrix D:\n", D)
    print(f"Symmetric: {is_symmetric(D)}")
    # The determinant is pi^3 - 1 > 0, and diagonal elements are positive. So it is PSD.
    print(f"PSD: True (det=pi^3-1 > 0)")
    print(f"Trace: {np.trace(D)} >= 1")
    print("D satisfies the necessary conditions.")
    print("However, D's entries involve the transcendental number pi.")
    print("If D were in P, it could be written as a finite convex combination of matrices with integer entries.")
    print("This would imply an algebraic relation for pi, e.g., of the form c2*pi^2 + c1*pi + c0 = 0 for integers c0, c1, c2, not all zero.")
    print("Since pi is transcendental, no such relation exists.")
    print("Conclusion: D is not in P.\n")

    # Matrix E
    print("\nMatrix E:\n", E)
    print(f"Symmetric: {is_symmetric(E)}")
    det_E = 1 - pi**2
    print(f"Determinant is 1 - pi^2 = {det_E} < 0. So E is not PSD.")
    print("Conclusion: E is not in P.\n")

    # Matrix F
    print("\nMatrix F:\n", F)
    print(f"Symmetric: {is_symmetric(F)}")
    print(f"PSD: {is_psd(F)}")
    print(f"Trace: {np.trace(F)} >= 1")
    print("F satisfies the necessary conditions. Let's try to construct it.")
    # Generators must have b=0, so they are of form [[a^2, 0], [0, 0]].
    # We need to express 42 as a convex combination of squares of integers.
    # Let's use 36 (6^2) and 49 (7^2).
    # 42 = l1*36 + l2*49. l1+l2=1.
    # l1*36 + (1-l1)*49 = 42  => 49 - 13*l1 = 42 => 13*l1 = 7 => l1 = 7/13.
    # So l2 = 6/13.
    l1, l2 = 7/13, 6/13
    S1 = np.array([[36, 0], [0, 0]]) # from v=(6,0)
    S2 = np.array([[49, 0], [0, 0]]) # from v=(7,0)
    constructed_F = l1 * S1 + l2 * S2
    if np.allclose(F, constructed_F):
        print("F can be constructed as a convex combination:")
        print(f"F = {l1:.4f} * {S1.tolist()} + {l2:.4f} * {S2.tolist()}")
        print("The equation with fractions is:")
        print(f"F = (7/13) * {S1.tolist()} + (6/13) * {S2.tolist()}")
        print(f"(7/13)*36 + (6/13)*49 = ({7*36 + 6*49})/13 = ({252 + 294})/13 = 546/13 = 42.")
        print("Conclusion: F is in P.\n")
        in_P.append('F')
    else:
        print("Conclusion: F is not in P.\n")

    print("\n--- Final Answer ---")
    print(in_P)

if __name__ == '__main__':
    check_matrices()
    final_answer = ['C', 'F']
    print(f"\n<<<{final_answer}>>>")