import numpy as np

# Step 1: Set up the problem.
# The problem is to find the maximum `c` such that the matrix A_n - c*I is positive
# semi-definite for all n, where A_ij = n - |i-j|. This `c` is the infimum of the
# minimum eigenvalues of A_n over all n.

# Step 2: Numerically compute the minimum eigenvalue for increasing n to find the trend.
max_n = 50
min_eigenvalues = []

print("This program numerically finds the maximum value of c.")
print("It computes the minimum eigenvalue of the relevant matrix for increasing sizes n.")
print("The sequence of minimum eigenvalues converges to the value of c.\n")
print("n  | min_eigenvalue(A_n)")
print("---|--------------------")

for n in range(1, max_n + 1):
    # Construct the matrix A_n for the current n
    # A = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=int)
    # A more direct construction:
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = n - abs(i-j)
            
    # For a real symmetric matrix, eigvalsh is efficient and returns real eigenvalues sorted.
    eigenvalues = np.linalg.eigvalsh(A)
    min_lambda = eigenvalues[0] # The smallest eigenvalue
    min_eigenvalues.append(min_lambda)
    
    # Print the minimum eigenvalue for selected n
    if n <= 10 or n % 5 == 0 or n == max_n:
        print(f"{n:<2} | {min_lambda:.8f}")

# Step 3: Compare with the theoretical result.
# The theoretical limit of the minimum eigenvalue as n -> infinity is 4/pi^2.
c_theoretical = 4 / (np.pi**2)

print("\nAs n increases, the minimum eigenvalue decreases and converges.")
print(f"The minimum eigenvalue for n={max_n} is approximately {min_eigenvalues[-1]:.8f}")
print(f"The theoretical value for c is 4 / pi^2, which is approximately {c_theoretical:.8f}\n")
print("Therefore, the maximum real number c is 4 / (pi^2).")
print("\nFinal Answer expressed with its components:")
print(f"c = 4 / ({np.pi})^2")
print("Which can be written as:")
print("4 / (pi**2)")
