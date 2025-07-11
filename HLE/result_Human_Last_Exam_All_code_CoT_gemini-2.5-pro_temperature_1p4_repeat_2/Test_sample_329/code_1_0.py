import numpy as np
from fractions import Fraction

def analyze_matrices():
    """
    Analyzes which of the given matrices are in the set P.
    """
    final_list = []

    # Matrix A
    print("--- Analyzing Matrix A ---")
    A = np.array([[0, 0], [0, 0]])
    print("A =\n", A)
    trace_A = np.trace(A)
    print("Trace of A:", trace_A)
    if trace_A < 1:
        print("Result: A is not in P because its trace must be >= 1.\n")
    else:
        # This part will not be reached for A
        final_list.append('A')

    # Matrix B
    print("--- Analyzing Matrix B ---")
    B = np.array([[6, 4], [3, 7]])
    print("B =\n", B)
    if not np.all(B == B.T):
        print("Result: B is not in P because it is not symmetric.\n")
    else:
        final_list.append('B')

    # Matrix C
    print("--- Analyzing Matrix C ---")
    C = np.array([[1, -1/2], [-1/2, 1]])
    print("C =\n", C)
    # Check necessary conditions
    is_symmetric = np.all(C == C.T)
    is_psd = C[0,0] >= 0 and C[1,1] >= 0 and np.linalg.det(C) >= 0
    has_valid_trace = np.trace(C) >= 1
    if not (is_symmetric and is_psd and has_valid_trace):
         print("Result: C fails one of the necessary conditions.\n")
    else:
        print("C satisfies the necessary conditions (symmetric, PSD, trace >= 1).")
        print("Attempting to construct C from generators...")
        v1 = np.array([[1],[1]])
        X1 = v1 @ v1.T
        v2 = np.array([[1],[-1]])
        X2 = v2 @ v2.T
        # Solve C = l1*X1 + l2*X2 with l1+l2=1
        # From C[0,1] = l1*1 + l2*(-1) = -0.5
        # l1 - (1-l1) = -0.5  => 2*l1 - 1 = -0.5 => 2*l1 = 0.5 => l1 = 0.25
        l1 = Fraction(1, 4)
        l2 = 1 - l1
        constructed_C = float(l1) * X1 + float(l2) * X2
        if np.allclose(C, constructed_C) and l1 >= 0 and l2 >= 0:
            print(f"C can be constructed with lambda1={l1}, lambda2={l2}.")
            print("The generating matrices are from v1=(1,1) and v2=(1,-1).")
            print(f"Final equation: {l1} * [[{X1[0,0]},{X1[0,1]}],[{X1[1,0]},{X1[1,1]}]] + {l2} * [[{X2[0,0]},{X2[0,1]}],[{X2[1,0]},{X2[1,1]}]] = [[{C[0,0]},{C[0,1]}],[{C[1,0]},{C[1,1]}]]")
            print("Result: C is in P.\n")
            final_list.append('C')
        else:
            print("Result: C could not be constructed as a convex combination of simple generators.\n")

    # Matrix D
    print("--- Analyzing Matrix D ---")
    pi = np.pi
    D = np.array([[pi, 1], [1, pi**2]])
    print("D =\n", D)
    print("A matrix in P lies in the affine hull of its integer-based generators.")
    print("This hull is a rational affine subspace (defined by linear equations with rational coefficients).")
    print("The point for D is (pi, 1, pi^2). A point with transcendental coordinates can't lie in a proper rational affine subspace.")
    print("For it to lie in such a subspace, a polynomial equation in pi with rational coefficients would have to be zero.")
    print("This contradicts the fact that pi is a transcendental number.")
    print("Result: D is not in P.\n")


    # Matrix E
    print("--- Analyzing Matrix E ---")
    E = np.array([[1, pi], [pi, 1]])
    print("E =\n", E)
    det_E = np.linalg.det(E)
    print(f"Determinant of E is 1 - pi^2 = {det_E:.4f}")
    if det_E < 0:
        print("Result: E is not in P because it is not positive semidefinite.\n")
    else:
        final_list.append('E')


    # Matrix F
    print("--- Analyzing Matrix F ---")
    F = np.array([[42, 0], [0, 0]])
    print("F =\n", F)
    is_symmetric = np.all(F == F.T)
    is_psd = F[0,0] >= 0 and F[1,1] >= 0 and np.linalg.det(F) >= 0
    has_valid_trace = np.trace(F) >= 1
    if not (is_symmetric and is_psd and has_valid_trace):
         print("Result: F fails one of the necessary conditions.\n")
    else:
        print("F satisfies the necessary conditions.")
        print("The (2,2) entry is 0, which implies only generators from v=(k,0) can be used.")
        print("We need to show 42 is a convex combination of integer squares k^2.")
        print("We can use k1^2 = 36 (from v=(6,0)) and k2^2 = 49 (from v=(7,0)).")
        # Solve 42 = l1*36 + l2*49 with l1+l2=1
        # 42 = l1*36 + (1-l1)*49 = 49 - 13*l1 => 13*l1 = 7
        l1 = Fraction(7, 13)
        l2 = 1 - l1
        v1 = np.array([[6],[0]])
        X1 = v1 @ v1.T
        v2 = np.array([[7],[0]])
        X2 = v2 @ v2.T
        constructed_F = float(l1) * X1 + float(l2) * X2
        if np.allclose(F, constructed_F) and l1 >= 0 and l2 >= 0:
            print(f"F can be constructed with lambda1={l1}, lambda2={l2}.")
            print(f"Final equation: {l1} * [[{X1[0,0]},{X1[0,1]}],[{X1[1,0]},{X1[1,1]}]] + {l2} * [[{X2[0,0]},{X2[0,1]}],[{X2[1,0]},{X2[1,1]}]] = [[{F[0,0]},{F[0,1]}],[{F[1,0]},{F[1,1]}]]")
            print("Result: F is in P.\n")
            final_list.append('F')
        else:
            print("Result: F could not be constructed.\n")

    print("--- Summary ---")
    print(f"The list of matrices contained in P is: {final_list}")
    
    # Final answer block
    print(f"<<<[{','.join(final_list)} T_T]>>>".replace(" ","").replace("T_T",""))


analyze_matrices()