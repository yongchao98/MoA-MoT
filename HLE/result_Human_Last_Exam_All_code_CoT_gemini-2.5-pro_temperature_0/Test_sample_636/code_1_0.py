import numpy as np

# --- Introduction ---
# We are given the system matrices A and B.
# A = [[-1, 1], [1, 0]]
# B = [[1, 2], [1, 0]]
# We want to find a feedback gain matrix F such that the eigenvalues of the
# closed-loop system A + B*F are at -1 + j and -1 - j.

# --- Derivation Summary ---
# 1. The desired characteristic polynomial for eigenvalues -1 ± j is s^2 + 2s + 2 = 0.
# 2. The characteristic polynomial of A + B*F is found by calculating det(sI - (A + B*F)).
#    Let F = [[f11, f12], [f21, f22]].
#    This yields: s^2 + s*(1 - f11 - f12 - 2*f21) + (-1 - f11 - 2*f12 + 2*f12*f21 - 2*f22 - 2*f11*f22) = 0.
# 3. Equating the coefficients with the desired polynomial gives two equations:
#    (I)  f11 + f12 + 2*f21 = -1
#    (II) -f11 - 2*f12 + 2*f12*f21 - 2*f22 - 2*f11*f22 = 3
# 4. To find a simple solution, we set f21 = 0 and f22 = 0. The system becomes:
#    (I)  f11 + f12 = -1
#    (II) -f11 - 2*f12 = 3
# 5. Solving this linear system gives f11 = 1 and f12 = -2.
# 6. Thus, a valid feedback gain matrix is F = [[1, -2], [0, 0]].

# --- Calculation and Verification ---

# Define system matrices
A = np.array([[-1., 1.],
              [1., 0.]])
B = np.array([[1., 2.],
              [1., 0.]])

# The calculated feedback gain matrix from our derivation
F = np.array([[1., -2.],
              [0.,  0.]])

# Calculate the closed-loop system matrix A_cl = A + B*F
A_cl = A + B @ F

# --- Output ---

print("The final equation A + B*F = A_cl is shown below with the calculated matrices.")

print("\nA =")
print(A)

print("\nB =")
print(B)

print("\nCalculated F =")
print(F)

print("\nResulting A + B*F =")
print(A_cl)

# Verification step: Calculate eigenvalues of the resulting closed-loop system
eigenvalues, _ = np.linalg.eig(A_cl)

print("\n--- Verification ---")
print("The eigenvalues of the resulting matrix A + B*F are:")
# Use np.round to handle potential floating point inaccuracies
print(np.round(eigenvalues, 5))
print("These match the desired eigenvalues of -1 ± j.")
