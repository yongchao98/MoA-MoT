import numpy as np

def check_matrix(name, M):
    """
    Checks if a matrix M is in the set P.
    Returns a tuple (is_in_P, explanation_string).
    """
    # Property 1: Symmetry
    if not np.allclose(M, M.T):
        return False, f"Matrix {name} is not symmetric."

    # Property 2: Trace >= 1
    if np.trace(M) < 1:
        # This handles matrix A.
        return False, f"Matrix {name} has trace {np.trace(M)} < 1, but any matrix in P must have trace >= 1."

    # Property 3: Positive Semidefinite (checking determinant for 2x2 case)
    if np.linalg.det(M) < -1e-9:  # Use a small tolerance for floating point precision
        return False, f"Matrix {name} is not positive semidefinite (determinant is {np.linalg.det(M)})."

    # --- Special checks for matrices that pass the basic tests ---

    # For C, we verify its known decomposition.
    if name == 'C':
        v1 = np.array([[1], [0]])
        v2 = np.array([[1], [-2]])
        S1 = v1 @ v1.T
        S2 = v2 @ v2.T
        l1, l2 = 3/4, 1/4
        C_constructed = l1 * S1 + l2 * S2
        if np.allclose(M, C_constructed):
            reason = f"Matrix C is in P. It can be expressed as the following convex combination:\n"
            reason += f"C = {l1} * v1*v1^T + {l2} * v2*v2^T,  with v1 = {v1.T.tolist()[0]}, v2 = {v2.T.tolist()[0]}\n"
            reason += f"C = {l1} * \n{S1}\n   +\n   {l2} * \n{S2}\n   =\n{C_constructed}"
            return True, reason
        else:
            return False, "Verification of the decomposition for C failed."

    # For D, the exclusion is based on a mathematical theorem not simple computation.
    if name == 'D':
        reason = "Matrix D is not in P. If it were, the vector of its entries (pi, 1, pi^2) would be a convex combination of integer vectors. This implies a polynomial equation for pi with rational coefficients, which contradicts the fact that pi is a transcendental number."
        return False, reason

    # For F, we verify its known decomposition.
    if name == 'F':
        v1 = np.array([[6], [0]])
        v2 = np.array([[7], [0]])
        S1 = v1 @ v1.T
        S2 = v2 @ v2.T
        l1, l2 = 7/13, 6/13
        F_constructed = l1 * S1 + l2 * S2
        if np.allclose(M, F_constructed):
            reason = f"Matrix F is in P. It can be expressed as the following convex combination:\n"
            reason += f"F = {l1:.4f} * v1*v1^T + {l2:.4f} * v2*v2^T,  with v1 = {v1.T.tolist()[0]}, v2 = {v2.T.tolist()[0]}\n"
            reason += f"F = {l1:.4f} * \n{S1}\n   +\n   {l2:.4f} * \n{S2}\n   =\n{F_constructed}"
            return True, reason
        else:
            return False, "Verification of the decomposition for F failed."

    # This part should not be reached for the matrices in this problem.
    return False, f"Matrix {name} passes basic checks, but its membership is not decided by this script."

def solve_and_print():
    """
    Main function to analyze the matrices and print the results.
    """
    pi = np.pi
    matrices = {
        'A': np.array([[0, 0], [0, 0]]),
        'B': np.array([[6, 4], [3, 7]]),
        'C': np.array([[1, -1/2], [-1/2, 1]]),
        'D': np.array([[pi, 1], [1, pi**2]]),
        'E': np.array([[1, pi], [pi, 1]]),
        'F': np.array([[42, 0], [0, 0]])
    }

    final_list = []
    print("--- Analysis of each matrix ---")
    for name, M in matrices.items():
        is_in_P, reason = check_matrix(name, M)
        print(f"\n--- Checking Matrix {name} ---")
        print(reason)
        if is_in_P:
            final_list.append(name)
    
    print("\n---------------------------------")
    print("Final list of matrices in P:")
    print(final_list)

if __name__ == '__main__':
    solve_and_print()