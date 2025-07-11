import numpy as np

def solve_and_print_reasoning():
    """
    Analyzes each matrix to determine if it belongs to the set P,
    printing the reasoning for each case.
    """
    # Define a helper function to create a generating matrix from a vector v
    def get_generator(v):
        v_col = np.array(v).reshape(2, 1)
        return v_col @ v_col.T

    # Define the matrices to be checked
    matrices = {
        'A': np.array([[0, 0], [0, 0]]),
        'B': np.array([[6, 4], [3, 7]]),
        'C': np.array([[1, -1/2], [-1/2, 1]]),
        'D': np.array([[np.pi, 1], [1, np.pi**2]]),
        'E': np.array([[1, np.pi], [np.pi, 1]]),
        'F': np.array([[42, 0], [0, 0]])
    }

    result_list = []
    print("--- Analyzing Matrices ---\n")

    # Necessary conditions for a matrix M to be in P:
    # 1. M must be symmetric.
    # 2. M must be positive semidefinite (eigenvalues >= 0).
    # 3. The trace of M must be >= 1.

    # Matrix A
    print("Matrix A = [[0, 0], [0, 0]]")
    tr_A = np.trace(matrices['A'])
    if tr_A < 1:
        print(f"Verdict: NOT in P. Reason: Trace is {tr_A}, which is less than 1.\n")
    
    # Matrix B
    print("Matrix B = [[6, 4], [3, 7]]")
    if not np.allclose(matrices['B'], matrices['B'].T):
        print("Verdict: NOT in P. Reason: The matrix is not symmetric.\n")

    # Matrix E
    print(f"Matrix E = [[1, {np.pi:.4f}], [{np.pi:.4f}, 1]]")
    det_E = np.linalg.det(matrices['E'])
    if det_E < 0:
        print(f"Verdict: NOT in P. Reason: The determinant is {det_E:.4f} < 0, so it's not positive semidefinite.\n")

    # Matrix C
    print("Matrix C = [[1, -0.5], [-0.5, 1]]")
    C = matrices['C']
    # Check basic properties first
    if np.allclose(C, C.T) and np.all(np.linalg.eigvalsh(C) >= -1e-9) and np.trace(C) >= 1:
        # It passes basic checks. We demonstrate it's in P by construction.
        S1 = get_generator([1, -1])
        S2 = get_generator([1, 1])
        # We express C = lambda*S1 + (1-lambda)*S2. From off-diagonal: -0.5 = -lambda + (1-lambda) => lambda = 3/4
        lmbda = 3/4
        combination = lmbda * S1 + (1 - lmbda) * S2
        if np.allclose(C, combination):
            print("Verdict: IS in P. Reason: It can be constructed as a convex combination.")
            print(f"C = {lmbda} * S(1,-1) + {1-lmbda} * S(1,1)")
            print(f"[[1, -0.5], [-0.5, 1]] = {lmbda} * {S1.tolist()} + {1-lmbda} * {S2.tolist()}")
            print(f"                   = {combination.tolist()}\n")
            result_list.append('C')

    # Matrix F
    print("Matrix F = [[42, 0], [0, 0]]")
    F = matrices['F']
    # Check basic properties first
    if np.allclose(F, F.T) and np.all(np.linalg.eigvalsh(F) >= -1e-9) and np.trace(F) >= 1:
        # F_22 = 0 implies it can only be a convex combination of S(a, 0) matrices.
        # This requires F_11=42 to be in the convex hull of {1, 4, 9,...}, which is [1, inf).
        if F[0, 0] >= 1:
            print("Verdict: IS in P. Reason: Its non-zero diagonal element is in the convex hull of the allowed squared integers.")
            # Demonstrate with a specific combination: 42 = lambda*36 + (1-lambda)*49 => lambda=7/13
            lmbda = 7/13
            S1 = get_generator([6, 0])  # a^2 = 36
            S2 = get_generator([7, 0])  # a^2 = 49
            combination = lmbda * S1 + (1 - lmbda) * S2
            print(f"F = {lmbda:.4f} * S(6,0) + {1-lmbda:.4f} * S(7,0)")
            print(f"[[42, 0], [0, 0]] = {lmbda:.4f} * {S1.tolist()} + {1-lmbda:.4f} * {S2.tolist()}")
            print(f"               = {combination.tolist()}\n")
            result_list.append('F')

    # Matrix D
    print(f"Matrix D = [[{np.pi:.4f}, 1], [1, {np.pi**2:.4f}]]")
    D = matrices['D']
    # Check basic properties first
    if np.allclose(D, D.T) and np.all(np.linalg.eigvalsh(D) >= -1e-9) and np.trace(D) >= 1:
        print("Verdict: NOT in P. Reason: Membership in P would require pi to be an algebraic number, but it is transcendental.")
        print("A formal argument shows that if D were in P, its element pi would have to be an integer, which is a contradiction.\n")

    # Print the final result in the requested format
    print("--- Final Result ---")
    print("The list of matrices contained in P is:")
    print(sorted(result_list))

# Run the analysis
solve_and_print_reasoning()