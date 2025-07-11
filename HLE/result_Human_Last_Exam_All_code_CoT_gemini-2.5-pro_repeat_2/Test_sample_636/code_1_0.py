import numpy as np
import control as ct

def solve_pole_placement():
    """
    Calculates the state feedback gain F for a given system so that
    the closed-loop system A + BF has specified eigenvalues.
    """
    # Define the system matrices A and B
    A = np.array([[-1., 1.],
                  [1., 0.]])
    B = np.array([[1., 2.],
                  [1., 0.]])

    # Define the desired eigenvalues (poles) for the closed-loop system
    desired_poles = np.array([-1. + 1.j, -1. - 1.j])

    # The control.place function solves for K in the equation A - BK.
    # Our problem is formulated as A + BF.
    # Therefore, we can find K and then F = -K.
    try:
        K = ct.place(A, B, desired_poles)
        F = -K
    except Exception as e:
        print(f"Error during pole placement: {e}")
        print("This may happen if the system is not controllable.")
        return

    # Calculate the resulting closed-loop matrix A_cl = A + B * F
    A_cl = A + B @ F

    # --- Output the results ---
    print("Given system matrix A:")
    print(A)

    print("\nGiven input matrix B:")
    print(B)

    print("\nDesired eigenvalues:")
    print(desired_poles)

    print("\nCalculated state feedback gain matrix F:")
    print(F)

    print("\nThe final equation is A_cl = A + B * F, where A_cl is the closed-loop matrix.")
    print("\nResulting closed-loop matrix A_cl:")
    print(A_cl)

    # Verification step: Calculate the eigenvalues of the new system A + BF
    eigenvalues_of_A_cl = np.linalg.eigvals(A_cl)
    print("\nVerification: Eigenvalues of the closed-loop system A_cl:")
    # Sort for consistent comparison
    eigenvalues_of_A_cl.sort(key=np.imag)
    print(eigenvalues_of_A_cl)

if __name__ == '__main__':
    solve_pole_placement()