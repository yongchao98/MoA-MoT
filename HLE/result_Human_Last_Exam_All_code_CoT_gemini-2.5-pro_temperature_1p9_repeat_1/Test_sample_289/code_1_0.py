import numpy as np

# This script verifies the result by constructing a matrix A that has the
# maximum possible number of non-real eigenvalues.

# From the derivation, any eigenvalue 'lambda' of a matrix A where A^3 = A*
# must satisfy the equation lambda^3 = conj(lambda).
# The solutions to this equation are lambda in {0, 1, -1, i, -i}.
# The set of possible non-real eigenvalues is therefore {i, -i}.
# This means the maximum size |S| is at most 2.

# Let's show this maximum is achievable by constructing a matrix A whose
# eigenvalues are exactly {i, -i}. The simplest choice is a diagonal matrix.
eigenvalues_to_test = np.array([1j, -1j], dtype=complex)
A = np.diag(eigenvalues_to_test)

# Verify the condition A^3 = A* for our constructed matrix A.
A_cubed = np.linalg.matrix_power(A, 3)
A_star = A.conj().T

# Check if the condition holds true
condition_holds = np.allclose(A_cubed, A_star)

# Find the set S of non-real eigenvalues of our matrix A
eigenvalues_of_A = np.linalg.eigvals(A)
S = {val for val in eigenvalues_of_A if not np.isclose(val.imag, 0)}

# Print the final result and the verification steps.
print("Step 1: The equation for the eigenvalues")
print("Any eigenvalue lambda must satisfy the equation: lambda**3 = conjugate(lambda)")
print(f"The possible non-real solutions are: {{1j, -1j}}")
print("This implies the largest possible size |S| is 2.")
print("-" * 30)

print("Step 2: Verification by construction")
print(f"We construct a matrix A = diag{tuple(eigenvalues_to_test)}.")
print(f"A =\n{A}")
print(f"We check if A**3 = A*:")
print(f"A**3 =\n{A_cubed}")
print(f"A* =\n{A_star}")
print(f"Is A**3 == A*? {condition_holds}")
print("-" * 30)

print("Step 3: Conclusion")
print(f"The set of non-real eigenvalues for our matrix A is S = {S}.")
print(f"The size of this set is |S| = {len(S)}.")
print("Since the maximum possible size is 2 and we have found a matrix that achieves it,")
print("the largest size |S| is 2.")
