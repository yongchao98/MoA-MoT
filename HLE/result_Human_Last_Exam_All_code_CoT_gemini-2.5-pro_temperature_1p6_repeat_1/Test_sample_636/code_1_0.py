import numpy as np

def solve_feedback_gain():
    """
    Calculates the state feedback gain F for a given linear system A, B
    such that the closed-loop system A + BF has specified eigenvalues.
    """
    # System matrices
    A = np.array([[-1, 1], 
                  [1, 0]])
    B = np.array([[1, 2], 
                  [1, 0]])

    # Desired eigenvalues (poles)
    desired_poles = np.array([-1 + 1j, -1 - 1j])

    # The desired characteristic polynomial is (lambda - p1)(lambda - p2)
    # For p1, p2 = -1 +/- j, this is lambda^2 + 2*lambda + 2 = 0
    #
    # The characteristic polynomial of A_cl = A + BF is det(lambda*I - A_cl) = 0
    # Let F = [[f1, f2], [f3, f4]].
    # A_cl = A + BF = [[-1+f1+2*f3, 1+f2+2*f4], [1+f1, f2]]
    # The polynomial is lambda^2 - trace(A_cl)*lambda + det(A_cl) = 0.
    #
    # Equating coefficients with the desired polynomial gives:
    # 1) -trace(A_cl) = 2  => -(-1+f1+2*f3 + f2) = 2  => f1+f2+2*f3 = -1
    # 2) det(A_cl) = 2     => f2*(-1+f1+2*f3) - (1+f1)*(1+f2+2*f4) = 2 => -f1-2*f2+2*f2*f3-2*f4-2*f1*f4 = 3
    #
    # This system is underdetermined. We find a simple solution by setting f3=0, f4=0.
    # The equations simplify to:
    # 1) f1 + f2 = -1
    # 2) -f1 - 2*f2 = 3
    #
    # This can be solved as a linear system M*x = v where x = [f1, f2].
    
    M = np.array([[1, 1], 
                  [-1, -2]])
    v = np.array([-1, 3])

    # Solve for f1 and f2
    try:
        f1_f2 = np.linalg.solve(M, v)
        f1 = f1_f2[0]
        f2 = f1_f2[1]
        
        # We assumed f3 and f4 are zero for this solution
        f3 = 0.0
        f4 = 0.0

        # Construct the final gain matrix F
        F = np.array([[f1, f2], [f3, f4]])

        print("The system matrices are:")
        print("A =")
        print(A)
        print("\nB =")
        print(B)

        print("\nThe calculated state feedback gain matrix F is:")
        print(F)

        # --- Verification ---
        A_cl = A + B @ F
        eigenvalues = np.linalg.eigvals(A_cl)

        print("\n--- Verification ---")
        print("The final equation for the closed-loop matrix is A_cl = A + B@F")
        print("\nA = \n{}\n".format(A))
        print("B = \n{}\n".format(B))
        print("F = \n{}\n".format(F))
        print("Resulting A_cl = A + B@F = \n{}\n".format(A_cl))

        # Sort for consistent comparison
        eigenvalues_sorted = np.sort(eigenvalues)
        desired_poles_sorted = np.sort(desired_poles)

        print("The eigenvalues of the closed-loop system are: [{:.2f}, {:.2f}]".format(
            eigenvalues_sorted[0], eigenvalues_sorted[1]))
        print("The desired eigenvalues were: [{:.2f}, {:.2f}]".format(
            desired_poles_sorted[0], desired_poles_sorted[1]))

    except np.linalg.LinAlgError:
        print("Could not solve for the feedback gain with the chosen simplification.")

solve_feedback_gain()