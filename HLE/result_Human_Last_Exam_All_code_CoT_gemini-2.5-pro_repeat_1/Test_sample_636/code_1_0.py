import numpy as np

# Define the system matrices A and B
A = np.array([[-1, 1], 
              [1, 0]])
B = np.array([[1, 2], 
              [1, 0]])

# We need to solve a system of linear equations to find the elements of F.
# From the plan explained above, the equations are:
# f1 + f2 = -1
# f1 + 2*f2 = -3
# This can be written as M * f = v
M = np.array([[1, 1], 
              [1, 2]])
v = np.array([-1, -3])

# Solve for f = [f1, f2]
f_solution = np.linalg.solve(M, v)
f1 = f_solution[0]
f2 = f_solution[1]

# Construct the state feedback gain matrix F with the simplifying assumption F = [[f1, f2], [0, 0]]
F = np.array([[f1, f2], 
              [0, 0]])

# Calculate the closed-loop matrix A_cl = A + B @ F
A_cl = A + B @ F

# Calculate the eigenvalues of the closed-loop system to verify the result
eigenvalues, _ = np.linalg.eig(A_cl)

# --- Output the results ---
print("The calculated state feedback gain matrix F is:")
print(F)

print("\nThe final equation A + B @ F = A_cl is constructed as follows:")
print("\nA =")
print(A)
print("\nB =")
print(B)
print("\nF =")
print(F)

print("\nResulting closed-loop matrix A_cl = A + B @ F:")
print(A_cl)

print("\n--- Verification ---")
print("The eigenvalues of the closed-loop matrix A_cl are:")
# Using np.round to clean up floating point inaccuracies
print(np.round(eigenvalues, decimals=5))
print("These values match the desired poles at -1.0 + 1.0j and -1.0 - 1.0j.")
