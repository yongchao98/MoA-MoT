import numpy as np

# This script calculates and verifies the state feedback gain F for a given LTI system.
# The method used is pole placement, simplified for a MIMO system.

# Define the system matrices A and B.
A = np.array([[-1., 1.],
              [1., 0.]])

B = np.array([[1., 2.],
              [1., 0.]])

# Based on manual calculation, one possible feedback matrix F is:
F = np.array([[1., -2.],
              [0., 0.]])

# Calculate the closed-loop system matrix A_cl = A + B*F.
A_cl = A + B @ F

# Calculate the eigenvalues of the closed-loop system to verify the result.
eigenvalues, _ = np.linalg.eig(A_cl)

# Print the final results in a formatted equation.
print("The calculated state feedback gain matrix F is:")
print(F)
print("\nThe final equation A + B*F = A_cl is:\n")

# Print the equation A + B * F = A_cl with each number
print(f"   [ {A[0,0]:>5.2f}  {A[0,1]:>5.2f} ]   [ {B[0,0]:>5.2f}  {B[0,1]:>5.2f} ]   [ {F[0,0]:>5.2f}  {F[0,1]:>5.2f} ]   [ {A_cl[0,0]:>5.2f}  {A_cl[0,1]:>5.2f} ]")
print(f"   [ {A[1,0]:>5.2f}  {A[1,1]:>5.2f} ] + [ {B[1,0]:>5.2f}  {B[1,1]:>5.2f} ] * [ {F[1,0]:>5.2f}  {F[1,1]:>5.2f} ] = [ {A_cl[1,0]:>5.2f}  {A_cl[1,1]:>5.2f} ]")
print("\n")

print("Verification: Eigenvalues of the closed-loop system A + BF are:")
# Format the complex numbers for printing
# Using np.sort to ensure consistent order for comparison
sorted_eigenvalues = np.sort(eigenvalues)
formatted_eigenvalues = [f"{eig.real:.2f} + {eig.imag:.2f}j" if eig.imag >= 0 else f"{eig.real:.2f} - {-eig.imag:.2f}j" for eig in sorted_eigenvalues]
print(formatted_eigenvalues)
