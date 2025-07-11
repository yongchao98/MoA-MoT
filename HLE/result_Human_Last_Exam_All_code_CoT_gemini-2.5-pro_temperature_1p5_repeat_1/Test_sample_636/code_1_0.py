import numpy as np

def solve_state_feedback():
    """
    Calculates the state feedback gain F for a given system A, B
    and desired eigenvalues.
    """
    # System matrices
    A = np.array([[-1, 1],
                  [1, 0]])
    B = np.array([[1, 2],
                  [1, 0]])

    # The desired eigenvalues are -1 + j and -1 - j.
    # The corresponding characteristic polynomial is s^2 + 2s + 2 = 0.

    # We assume a simplified F of the form [[f1, f2], [0, 0]].
    # Comparing coefficients of the characteristic polynomial of A+BF
    # with the desired one gives a system of linear equations:
    # 1 - f1 - f2 = 2  =>  f1 + f2 = -1
    # -1 - f1 - 2*f2 = 2  =>  f1 + 2*f2 = -3

    # We solve this system M * x = y for x = [f1, f2]
    M = np.array([[1, 1],
                  [1, 2]])
    y = np.array([-1, -3])

    # Solve for f1 and f2
    f1_f2 = np.linalg.solve(M, y)
    f1 = f1_f2[0]
    f2 = f1_f2[1]

    # Construct the full gain matrix F
    F = np.array([[f1, f2],
                  [0, 0]])

    print("Calculated state feedback gain matrix F:")
    print(F)

    # Verification step
    A_closed_loop = A + B @ F
    eigenvalues = np.linalg.eigvals(A_closed_loop)

    print("\nVerification:")
    print("The eigenvalues of the resulting matrix A + BF are:")
    # np.round is used to clean up floating point inaccuracies
    print(np.round(eigenvalues, 5))
    
    # Per the instruction to "output each number in the final equation",
    # we print the elements of the F matrix that forms our solution.
    print("\nThe final equation for the feedback gain is F = [[f11, f12], [f21, f22]], where:")
    print(f"f11 = {F[0, 0]}")
    print(f"f12 = {F[0, 1]}")
    print(f"f21 = {F[1, 0]}")
    print(f"f22 = {F[1, 1]}")

solve_state_feedback()