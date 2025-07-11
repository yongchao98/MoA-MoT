import numpy as np

# Set print options for better readability
np.set_printoptions(precision=4, suppress=True)

# Step 1: Define the system and input matrices A and B.
A = np.array([[-1, 1],
              [1, 0]])
B = np.array([[1, 2],
              [1, 0]])

# Step 2: Define the desired eigenvalues and find the target closed-loop matrix.
# The desired characteristic polynomial for eigenvalues -1 ± j is s^2 + 2s + 2 = 0.
# We choose a target matrix A_cl (in companion form) with this polynomial.
A_cl_desired = np.array([[0,   1],
                         [-2, -2]])

# Step 3: Solve the equation A + B*F = A_cl for F.
# This becomes F = inv(B) * (A_cl - A)
try:
    F = np.linalg.inv(B) @ (A_cl_desired - A)
except np.linalg.LinAlgError:
    print("Matrix B is not invertible. This specific method fails.")
    F = None

if F is not None:
    # Step 4: Print all the matrices and the final verification.
    print("--- System Matrices ---")
    print(f"A =\n{A}\n")
    print(f"B =\n{B}\n")
    
    print("--- Calculation ---")
    print(f"The calculated state feedback gain matrix F is:")
    print(F)
    print("\nThis F is the solution to the equation A + B*F = A_cl, where A_cl is:")
    print(A_cl_desired)

    print("\n--- Verification ---")
    # The resulting closed-loop system from A + B*F
    A_cl_final = A + B @ F

    # Calculate the eigenvalues of the final closed-loop system
    eigenvalues = np.linalg.eigvals(A_cl_final)
    
    print("The resulting closed-loop matrix A + B*F is:")
    print(A_cl_final)
    print("\nThe eigenvalues of the final closed-loop system A + B*F are:")
    print(eigenvalues)
    print("\nThese match the desired eigenvalues: [-1.+1.j, -1.-1.j]")
