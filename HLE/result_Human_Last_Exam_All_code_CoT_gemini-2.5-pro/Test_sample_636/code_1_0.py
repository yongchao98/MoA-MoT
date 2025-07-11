import numpy as np

def solve_pole_placement():
    """
    Calculates the state feedback gain F for a given system so that
    A + BF has specified eigenvalues.
    """
    # Step 1: Define the system matrices
    A = np.array([[-1, 1], [1, 0]])
    B = np.array([[1, 2], [1, 0]])

    # Step 2: Define desired eigenvalues and find the characteristic polynomial coefficients
    # Desired eigenvalues are -1 + j and -1 - j.
    # Desired characteristic polynomial is s^2 + 2s + 2 = 0.
    # This means the trace of the closed-loop matrix (A+BF) must be -2,
    # and its determinant must be 2.
    desired_trace = -2
    desired_det = 2

    # Step 3: Formulate equations for the elements of F.
    # Let F = [[f1, f2], [f3, f4]].
    # To find a unique solution, we can simplify F by setting some elements to zero.
    # Let's assume F = [[f1, f2], [0, 0]].
    #
    # The closed-loop matrix A_cl = A + B*F becomes:
    # A_cl = [[-1, 1], [1, 0]] + [[1, 2], [1, 0]] * [[f1, f2], [0, 0]]
    #      = [[-1, 1], [1, 0]] + [[f1, f2], [f1, f2]]
    #      = [[-1 + f1, 1 + f2], [1 + f1, f2]]
    #
    # Now, we use the trace and determinant conditions:
    # trace(A_cl) = (-1 + f1) + f2 = f1 + f2 - 1
    # det(A_cl) = (-1 + f1)*f2 - (1 + f2)*(1 + f1)
    #
    # This gives us a system of two linear equations:
    # 1) f1 + f2 - 1 = desired_trace  =>  f1 + f2 = -1
    # 2) (-1 + f1)*f2 - (1 + f2)*(1 + f1) = desired_det
    #    -f2 + f1*f2 - (1 + f1 + f2 + f1*f2) = 2
    #    -f2 + f1*f2 - 1 - f1 - f2 - f1*f2 = 2
    #    -f1 - 2*f2 - 1 = 2          => -f1 - 2*f2 = 3

    # We solve the system:
    #  1*f1 +  1*f2 = -1
    # -1*f1 + -2*f2 =  3
    
    # In matrix form M * x = v, where x = [f1, f2]
    M = np.array([[1, 1], [-1, -2]])
    v = np.array([-1, 3])

    # Step 4: Solve for f1 and f2
    try:
        f1_f2 = np.linalg.solve(M, v)
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # Construct the full F matrix with f3=0, f4=0
    F = np.zeros((2, 2))
    F[0, :] = f1_f2
    
    # Calculate the resulting closed-loop system for verification
    A_cl = A + B @ F

    # Step 5: Print the results and verification
    print("The calculated state feedback gain matrix F is:")
    print(F)
    print("\nTo verify the solution, we form the equation A + B*F = A_closed_loop:")
    
    # Format the matrices into strings for clear output
    A_str = np.array2string(A, precision=2, separator=', ')
    B_str = np.array2string(B, precision=2, separator=', ')
    F_str = np.array2string(F, precision=2, separator=', ')
    A_cl_str = np.array2string(A_cl, precision=2, separator=', ')
    
    print(f"{A_str}\n+ {B_str} * {F_str}\n= {A_cl_str}")

    # Final verification by calculating eigenvalues
    eigenvalues = np.linalg.eigvals(A_cl)
    print("\nThe eigenvalues of the resulting matrix A + BF are:")
    print(eigenvalues)

solve_pole_placement()