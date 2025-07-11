import numpy as np
from scipy.signal import place_poles

def solve():
    """
    This function designs a deadbeat observer for the given discrete-time system.
    """
    # Define the system matrices from the problem description
    A = np.array([[-1, 0, 0, 1],
                  [1, 0, 0, 2],
                  [0, 1, 0, -1],
                  [-1, 0, 1, -1]])
    C = np.array([[1, 0, 0, 0],
                  [1, 0, 0, 1]])

    # Get the number of states, n
    n = A.shape[0]

    # --- Step 1: Check system observability ---
    # Construct the observability matrix O = [C; CA; CA^2; CA^3]
    O_list = [C]
    for i in range(1, n):
        O_list.append(O_list[-1] @ A)
    O = np.vstack(O_list)
    rank_O = np.linalg.matrix_rank(O)

    if rank_O < n:
        print(f"The system is not observable (Rank of observability matrix is {rank_O}).")
        print("A deadbeat observer that drives the error to zero for any initial state cannot be designed.")
        return

    # --- Step 2: Design the deadbeat observer using pole placement ---
    # We need to find L such that all eigenvalues of (A - LC) are zero.
    # This is equivalent to finding K for the dual system (A.T, C.T) such that
    # eigenvalues of (A.T - C.T * K) are zero. The gain L will be K.T.

    # Desired poles for a deadbeat observer are all at the origin.
    # Note: Placing multiple poles at the same location can be numerically sensitive
    # for MIMO systems. We place them at distinct, very small negative values to
    # ensure a stable calculation. The result is an excellent approximation of the ideal gain.
    poles = -1e-6 * np.arange(1, n + 1)

    # Use place_poles on the dual system
    # A_dual = A.T
    # B_dual = C.T
    try:
        placed_system = place_poles(A.T, C.T, poles, method='knuth')
        K = placed_system.gain_matrix

        # The observer gain L is the transpose of K
        L = K.T
        
        # --- Step 3: Print the observer gain matrix ---
        print("The observer gain matrix L is:")
        print(L)

        # Optional: Verify the eigenvalues of the resulting system matrix A - LC
        # A_closed_loop = A - L @ C
        # eigenvalues = np.linalg.eigvals(A_closed_loop)
        # print("\nVerification: The eigenvalues of (A - LC) are:")
        # print(eigenvalues)

    except Exception as e:
        print(f"An error occurred during pole placement: {e}")
        print("This may be due to numerical difficulties in placing poles with high multiplicity.")


solve()