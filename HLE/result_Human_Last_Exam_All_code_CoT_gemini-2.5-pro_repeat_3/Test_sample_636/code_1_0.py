import numpy as np

def solve_and_verify_pole_placement():
    """
    Calculates the state feedback gain F for a given system
    such that A + BF has specified eigenvalues.
    """
    # Set numpy print options for clean output
    np.set_printoptions(precision=4, suppress=True)

    # Step 1: Define the system matrices A and B
    A = np.array([[-1., 1.],
                  [1., 0.]])
    B = np.array([[1., 2.],
                  [1., 0.]])

    # Step 2: Define the desired eigenvalues (poles)
    # The desired characteristic polynomial for poles at -1 +/- j is:
    # (s - (-1+j))(s - (-1-j)) = s^2 + 2s + 2 = 0
    poles = np.array([-1. + 1.j, -1. - 1.j])

    # Step 3: Calculate the state feedback gain matrix F.
    # For a multi-input system, the solution for F is not unique.
    # We find one solution by assuming F has a simpler form, e.g., F = [[f1, f2], [0, 0]].
    # By solving the characteristic polynomial equations, we get:
    # f1 + f2 = -1
    # f1 + 2*f2 = -3
    # Solving this system gives f1 = 1 and f2 = -2.
    F = np.array([[1., -2.],
                  [0., 0.]])

    print("The calculated state feedback gain matrix F is:")
    print(F)

    # Step 4: Form the closed-loop system matrix A_cl = A + B @ F
    # and print the full equation as requested.
    A_cl = A + B @ F

    print("\nThe final equation A + B * F = A_cl is:")
    print(f"\n{A}\n\n+ \n\n{B}\n\n@\n\n{F}\n\n= \n\n{A_cl}")

    # Step 5: Verify the eigenvalues of the closed-loop system
    eigenvalues, _ = np.linalg.eig(A_cl)

    print("\n----------------------------------------------------")
    print("Verification:")
    print("The eigenvalues of the resulting matrix A + B * F are:")
    # Sorting for consistent comparison
    print(np.sort(eigenvalues))
    print("\nThe desired eigenvalues were:")
    print(np.sort(poles))
    print("----------------------------------------------------")

    # Check if the real and imaginary parts of the eigenvalues match
    if np.allclose(np.sort(eigenvalues), np.sort(poles)):
        print("\nVerification successful: The eigenvalues match the desired poles.")
    else:
        print("\nVerification failed: The eigenvalues do not match the desired poles.")

# Run the function
solve_and_verify_pole_placement()