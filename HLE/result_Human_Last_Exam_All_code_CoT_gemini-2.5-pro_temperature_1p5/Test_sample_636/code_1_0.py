import numpy as np

# Set print options for better readability
np.set_printoptions(precision=4, suppress=True)

# 1. Define the system matrices and desired eigenvalues
A = np.array([[-1., 1.], 
              [1., 0.]])
B = np.array([[1., 2.], 
              [1., 0.]])
desired_poles = np.array([-1. + 1j, -1. - 1j])

print("System Matrices:")
print("A = \n", A)
print("B = \n", B)
print("\nDesired eigenvalues (poles): ", desired_poles)

# 2. Determine the desired characteristic polynomial: s^2 + 2s + 2 = 0
# The coefficients from s^2 - trace(A+BF)s + det(A+BF) = 0 must match.
# -trace(A+BF) = 2  => 1 - f1 - f2 - 2*f3 = 2 => f1 + f2 + 2*f3 = -1
# det(A+BF) = 2 => -1 - f1 - 2*f2 + 2*f2*f3 - 2*f4(1+f1) = 2 => -f1 - 2*f2 + 2*f2*f3 - 2*f4 - 2*f1*f4 = 3

# 3. Simplify by setting f3=0 and f4=0. The system of equations becomes:
# f1 + f2 = -1
# -f1 - 2*f2 = 3
# We can write this as M * x = v, where x = [f1, f2]
M = np.array([[1., 1.], 
              [-1., -2.]])
v = np.array([-1., 3.])

# 4. Solve for f1 and f2
try:
    f1_f2 = np.linalg.solve(M, v)
    f1, f2 = f1_f2
    f3, f4 = 0., 0.

    # 5. Construct the feedback gain matrix F
    F = np.array([[f1, f2], [f3, f4]])

    print("\nCalculated State Feedback Gain Matrix F:")
    print("F = \n", F)

    # 6. Verification
    # Calculate the closed-loop system matrix A_cl = A + BF
    A_cl = A + B @ F

    print("\nResulting Closed-Loop Matrix (A + BF):")
    # This loop prints each number in the final equation matrix
    for row in A_cl:
        for val in row:
            print(f"{val:8.4f}", end="")
        print()


    # Calculate the eigenvalues of the closed-loop system
    eigenvalues = np.linalg.eigvals(A_cl)

    print("\nVerification: Eigenvalues of (A + BF):")
    print(eigenvalues)

    # Check if the calculated eigenvalues match the desired ones
    if np.allclose(sorted(eigenvalues.real), sorted(desired_poles.real)) and \
       np.allclose(sorted(abs(eigenvalues.imag)), sorted(abs(desired_poles.imag))):
        print("\nVerification successful: The eigenvalues match the desired poles.")
    else:
        print("\nVerification failed: The eigenvalues do not match the desired poles.")

except np.linalg.LinAlgError:
    print("\nThe system of equations for F has no unique solution or is singular.")
