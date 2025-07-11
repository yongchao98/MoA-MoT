import numpy as np

def solve_state_feedback():
    """
    Calculates the state feedback gain F for a given system (A, B)
    to place the closed-loop eigenvalues at desired locations.
    """
    # Define the system matrices
    A = np.array([[-1, 1], [1, 0]])
    B = np.array([[1, 2], [1, 0]])

    # Define the desired eigenvalues (poles)
    # The desired characteristic polynomial is s^2 + 2s + 2 = 0.
    
    # We need to find F = [[f1, f2], [f3, f4]] such that A+BF has the desired eigenvalues.
    # This leads to 2 equations for 4 unknowns, so there are infinite solutions.
    # We will find one solution by assuming f3=0 and f4=0.
    
    # The characteristic polynomial of A+BF with F = [[f1, f2], [0, 0]] is:
    # s^2 - (f1 + f2 - 1)s + (-1 - f1 - 2*f2) = 0
    
    # By comparing coefficients with the desired s^2 + 2s + 2 = 0, we get:
    # 1) -(f1 + f2 - 1) = 2  =>  f1 + f2 = -1
    # 2) -1 - f1 - 2*f2 = 2  =>  f1 + 2*f2 = -3
    
    # We solve this system of linear equations: M * x = v
    M = np.array([[1, 1], [1, 2]])
    v = np.array([-1, -3])
    
    # Solve for x = [f1, f2]
    f1_f2 = np.linalg.solve(M, v)
    
    # Construct the full gain matrix F
    F = np.array([[f1_f2[0], f1_f2[1]], [0.0, 0.0]])
    
    print("The calculated state feedback gain matrix F is:")
    print(F)
    
    # Verification
    A_cl = A + B @ F
    eigenvalues = np.linalg.eigvals(A_cl)
    
    print("\n--- Verification ---")
    print("The final equation A + B*F = A_cl is:")
    print(f"{A}\n+\n{B} @ {F}\n=\n{A_cl}")
    
    print("\nThe eigenvalues of the resulting closed-loop system A + BF are:")
    # Use np.round to clean up floating point inaccuracies
    print(np.round(eigenvalues, 5))
    print("These match the desired eigenvalues of -1+j and -1-j.")

solve_state_feedback()
<<<[[1. -2.]
 [0.  0.]]>>>