import numpy as np

def solve():
    """
    Determines which of the given matrices are in the set P,
    where P = ConvexHull(v * v^T : v in Z^2 \\ {(0,0)}).
    The final list of matrices in P is printed.
    """
    
    # Define the matrices using floating point numbers for calculations.
    pi = np.pi
    A = np.array([[0, 0], [0, 0]], dtype=float)
    B = np.array([[6, 4], [3, 7]], dtype=float)
    C = np.array([[1, -0.5], [-0.5, 1]], dtype=float)
    D = np.array([[pi, 1], [1, pi**2]])
    E = np.array([[1, pi], [pi, 1]])
    F = np.array([[42, 0], [0, 0]], dtype=float)

    matrices = {'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F}
    result_list = []
    
    print("--- Analysis of Matrices ---")

    # Check A
    if np.trace(A) < 1:
        print("A is not in P (trace is not >= 1).")
    
    # Check B
    if not np.allclose(B, B.T):
        print("B is not in P (not symmetric).")

    # Check E
    if np.linalg.det(E) < 0:
        print("E is not in P (not positive semidefinite, determinant is negative).")

    # Check C
    # C = (3/4) * M_v1 + (1/4) * M_v2 for v1=(1,0), v2=(1,-2)
    lambda1_C, lambda2_C = 3/4, 1/4
    v1_C = np.array([[1], [0]])
    v2_C = np.array([[1], [-2]])
    constructed_C = lambda1_C * (v1_C @ v1_C.T) + lambda2_C * (v2_C @ v2_C.T)
    if np.allclose(C, constructed_C):
        print("C is in P.")
        print(f"Decomposition: C = {lambda1_C}*[[1,0],[0,0]] + {lambda2_C}*[[1,-2],[-2,4]]")
        result_list.append('C')

    # Check D
    print("D is not in P (due to transcendental entries, as explained in the text).")

    # Check F
    # F = (7/13) * M_v1 + (6/13) * M_v2 for v1=(6,0), v2=(7,0)
    lambda1_F, lambda2_F = 7/13, 6/13
    v1_F = np.array([[6], [0]])
    v2_F = np.array([[7], [0]])
    constructed_F = lambda1_F * (v1_F @ v1_F.T) + lambda2_F * (v2_F @ v2_F.T)
    if np.allclose(F, constructed_F):
        print("F is in P.")
        print(f"Decomposition: F = (7/13)*[[36,0],[0,0]] + (6/13)*[[49,0],[0,0]]")
        print(f"Equation for F(1,1): (7/13)*36 + (6/13)*49 = {constructed_F[0,0]}")
        result_list.append('F')
        
    result_list.sort()
    
    print("\nFinal list of matrices in P:")
    print(result_list)

solve()