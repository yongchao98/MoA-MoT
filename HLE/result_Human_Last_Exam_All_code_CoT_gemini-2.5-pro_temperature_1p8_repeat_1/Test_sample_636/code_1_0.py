import numpy as np

def solve_and_verify():
    """
    Calculates the state feedback gain F for a given system to place
    eigenvalues at desired locations.
    """
    # Define system matrices
    A = np.array([[-1.0, 1.0], [1.0, 0.0]])
    B = np.array([[1.0, 2.0], [1.0, 0.0]])

    # The problem of placing eigenvalues at -1 +/- j leads to an
    # underdetermined system of equations for the elements of F.
    # By simplifying F to have a second row of zeros, we get a
    # solvable linear system for the first row [f1, f2]:
    #  1*f1 + 1*f2 = -1
    # -1*f1 - 2*f2 =  3

    # We can write this system in matrix form M*x = v, where x = [f1, f2]'
    M = np.array([[1.0, 1.0], [-1.0, -2.0]])
    v = np.array([-1.0, 3.0])

    # Solve for x = [f1, f2]
    try:
        solution = np.linalg.solve(M, v)
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # Construct the full gain matrix F from the solution
    F = np.array([
        [solution[0], solution[1]],
        [0.0, 0.0]
    ])

    # Calculate the closed-loop matrix A_cl = A + B*F for display
    A_cl = A + B @ F

    # Print all the matrices involved in the final equation as requested
    print("The goal is to find F such that the eigenvalues of A + B*F are at -1 +/- j.")
    print("A solution for F is:")
    print("F =\n", F)

    print("\nThe full equation A + B*F = A_cl contains the following matrices:")
    print("\nA =\n", A)
    print("\nB =\n", B)
    print("\nResulting in A_cl =\n", A_cl)


    # --- Verification Step ---
    eigenvalues = np.linalg.eigvals(A_cl)
    print("\n--- Verification ---")
    print("The eigenvalues of the final matrix A_cl are:")
    # Using np.round to handle potential floating point inaccuracies
    print(np.round(eigenvalues, 5))

# Execute the function
solve_and_verify()