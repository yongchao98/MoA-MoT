import numpy as np

# Set print options for better readability of matrices and complex numbers
np.set_printoptions(precision=3, suppress=True)

# 1. Define the system matrices A and B
A = np.array([
    [-1, 1],
    [1, 0]
])

B = np.array([
    [1, 2],
    [1, 0]
])

# 2. Define the calculated state feedback gain matrix F
# This F was calculated by hand to place the eigenvalues at -1 +/- j
F = np.array([
    [-3., 0.],
    [1., 0.]
])

# 3. Calculate the closed-loop system matrix A_cl = A + BF
A_cl = A + B @ F

print("Calculated Feedback Gain Matrix F:")
print(F)
print("\nResulting Closed-Loop Matrix A + BF:")
print(A_cl)

# 4. Verify the properties of the closed-loop matrix A_cl
# The characteristic equation is lambda^2 - trace(A_cl)*lambda + det(A_cl) = 0
# For desired eigenvalues -1 +/- j, the equation is lambda^2 + 2*lambda + 2 = 0
# This requires trace(A_cl) = -2 and det(A_cl) = 2.

# Calculate the trace and determinant of our resulting matrix
trace_A_cl = np.trace(A_cl)
det_A_cl = np.linalg.det(A_cl)

print("\n--- Verification of Characteristic Equation ---")
print("Desired characteristic equation: lambda^2 + 2*lambda + 2 = 0")
print(f"This requires trace = -2 and determinant = 2.")
print(f"Calculated trace(A+BF): {trace_A_cl:.1f}")
print(f"Calculated det(A+BF):   {det_A_cl:.1f}")

# 5. Calculate the eigenvalues of the closed-loop system to confirm the result
eigenvalues = np.linalg.eigvals(A_cl)

print("\nFinal Eigenvalues of A + BF:")
print(eigenvalues)
print("\nThe eigenvalues match the desired location of -1 +/- j.")
