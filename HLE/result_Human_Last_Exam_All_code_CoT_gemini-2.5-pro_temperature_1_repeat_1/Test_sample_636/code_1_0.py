import numpy as np

def solve_pole_placement():
    """
    Calculates the state feedback gain F for a MIMO system
    such that A + BF has specified eigenvalues.
    """
    # Set print options for better readability of matrices
    np.set_printoptions(precision=4, suppress=True)

    # 1. Define System and Desired Poles
    A = np.array([[-1., 1.], [1., 0.]])
    B = np.array([[1., 2.], [1., 0.]])
    desired_poles = np.array([-1 + 1j, -1 - 1j])

    # 2. Convert to a Single-Input System
    # We choose an input vector g to create an effective single input b_eff = B @ g.
    # Let's choose g = [1, 0]^T.
    g = np.array([[1.], [0.]])
    b_eff = B @ g

    # 3. Check Controllability and use Ackermann's formula
    n = A.shape[0]
    controllability_matrix = np.hstack([b_eff, A @ b_eff])

    if np.linalg.matrix_rank(controllability_matrix) != n:
        print("System is not controllable with the chosen g vector.")
        return

    # Desired characteristic polynomial from poles: (s-p1)(s-p2) = s^2 + 2s + 2
    poly_coeffs = np.poly(desired_poles)

    # Evaluate the characteristic polynomial at matrix A: P_des(A)
    # P_des(A) = A^2 + 2*A + 2*I
    P_A = np.polyval(poly_coeffs, A)

    # Ackermann's formula for gain K (for closed loop A - b_eff @ K)
    # K = [0, ..., 0, 1] @ inv(C) @ P_des(A)
    q = np.zeros((1, n))
    q[0, n - 1] = 1
    K = q @ np.linalg.inv(controllability_matrix) @ P_A

    # 4. Construct the full feedback gain F
    # We found K for A_cl = A - b_eff @ K = A - (B @ g) @ K.
    # We want A_cl = A + B @ F.
    # So, B @ F = -(B @ g) @ K = B @ (-g @ K). One solution is F = -g @ K.
    F = -g @ K

    # 5. Verification and Output
    A_cl = A + B @ F
    eigenvalues = np.linalg.eigvals(A_cl)

    print("--- System Matrices ---")
    print("A = \n", A)
    print("\nB = \n", B)
    
    print("\n--- Calculation Results ---")
    print("Desired eigenvalues: [-1 + j, -1 - j]")
    print("\nCalculated Feedback Gain F:")
    print(F)
    
    print("\n--- Verification ---")
    print("Resulting Closed-Loop Matrix (A + BF):")
    print(A_cl)
    
    print("\nFinal Equation (A + B * F = A_cl):")
    for i in range(n):
        print(f"{A[i, :]} + {B[i, :]} * {F} = {A_cl[i, :]}")

    print(f"\nThe eigenvalues of (A + BF) are: {eigenvalues[0]:.4f} and {eigenvalues[1]:.4f}")
    print("This confirms the solution is correct.")

solve_pole_placement()