import numpy as np

# Step 1: Define a matrix A whose eigenvalues are the set of non-real solutions {i, -i}
# We choose A = diag(i, -i).
A = np.array([[1j, 0], [0, -1j]])

print("Matrix A:")
print(A)

# Step 2: Calculate A^3 and A* to verify the condition A^3 = A*
A_cubed = np.linalg.matrix_power(A, 3)
A_star = A.conj().T

print("\nA^3:")
print(A_cubed)
print("\nA* (Adjoint of A):")
print(A_star)

# Check if A^3 is equal to A*
are_equal = np.allclose(A_cubed, A_star)
print(f"\nIs A^3 = A*? {are_equal}")

# Step 3: Find the eigenvalues of A
eigenvalues = np.linalg.eigvals(A)
print(f"\nEigenvalues of A: {eigenvalues}")

# Step 4: Identify the non-real eigenvalues and count them
non_real_eigenvalues = [val for val in eigenvalues if np.imag(val) != 0]
print(f"Non-real eigenvalues in S: {non_real_eigenvalues}")

# The size of the set S is the number of non-real eigenvalues
size_S = len(non_real_eigenvalues)
print(f"\nThe largest possible size |S| is {size_S}")
