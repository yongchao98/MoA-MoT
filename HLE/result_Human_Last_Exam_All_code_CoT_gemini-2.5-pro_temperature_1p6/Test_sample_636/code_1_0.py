import numpy as np

# Set print options for better readability
np.set_printoptions(precision=5, suppress=True)

# Define the system matrices A and B
A = np.array([
    [-1, 1],
    [1, 0]
])

B = np.array([
    [1, 2],
    [1, 0]
])

# Define the desired eigenvalues (poles)
desired_poles = np.array([-1 + 1j, -1 - 1j])

print("System Matrix A:")
print(A)
print("\nInput Matrix B:")
print(B)
print("\nDesired eigenvalues (poles):")
print(desired_poles)

# From our manual calculation, we found a possible F by assuming its second row is zero.
# f11 = 1, f12 = -2, f21 = 0, f22 = 0
F = np.array([
    [1.0, -2.0],
    [0.0, 0.0]
])

print("\nCalculated State Feedback Gain Matrix F:")
print(F)

# Calculate the closed-loop system matrix A_cl = A + BF
A_cl = A + B @ F

print("\nClosed-loop system matrix A + BF:")
print(A_cl)

# Calculate the eigenvalues of the closed-loop system to verify our solution
eigenvalues = np.linalg.eigvals(A_cl)

print("\nEigenvalues of the resulting closed-loop system A + BF:")
print(eigenvalues)

# Final verification
# The real parts should be close to -1 and imaginary parts close to +/- 1.
if np.allclose(np.sort(np.real(eigenvalues)), np.sort(np.real(desired_poles))) and \
   np.allclose(np.sort(np.abs(np.imag(eigenvalues))), np.sort(np.abs(np.imag(desired_poles)))):
    print("\nVerification successful: The calculated eigenvalues match the desired eigenvalues.")
else:
    print("\nVerification failed: The eigenvalues do not match.")
