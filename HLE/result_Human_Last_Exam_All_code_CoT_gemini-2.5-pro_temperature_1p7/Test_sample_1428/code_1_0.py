import numpy as np

# Set the number of sites on the circle
try:
    N_input = input("Enter the number of sites N (e.g., 10): ")
    N = int(N_input)
    if N < 3:
        print("N must be at least 3. Using N=10.")
        N = 10
except (ValueError, TypeError):
    print("Invalid input. Using default N=10.")
    N = 10

print("\n--- 1. One-dimensional Random Walk on a Circle ---")
print(f"Considering a circle with N = {N} sites.")
print("The probability distribution p(t) over the sites evolves as:")
print("p_i(t+1) = 0.5 * p_{i-1}(t) + 0.5 * p_{i+1}(t)\n")


print("--- 2. Transition Probability Matrix A ---")
print("This transformation is represented by the N x N matrix A:")
# Create the transition matrix A
A = np.zeros((N, N))
for i in range(N):
    # Probability 0.5 to move to the previous site (i-1)
    A[i, (i - 1 + N) % N] = 0.5
    # Probability 0.5 to move to the next site (i+1)
    A[i, (i + 1) % N] = 0.5

# Set numpy print options for better matrix display
if N > 10:
    np.set_printoptions(threshold=100, linewidth=120)
print(A)
print("\n")


print("--- 3. Eigenvalues of A ---")
print("The analytical eigenvalues are given by the formula: lambda_n = cos(2*pi*n/N)")
# Calculate analytical eigenvalues and sort them in descending order
analytical_eigenvalues = np.sort(np.cos(2 * np.pi * np.arange(N) / N))[::-1]

# Calculate numerical eigenvalues from the matrix A
numerical_eigenvalues, _ = np.linalg.eig(A)
# Sort the real parts in descending order for comparison
numerical_eigenvalues_sorted = np.sort(np.real(numerical_eigenvalues))[::-1]

print("\nNumerically computed eigenvalues of A (sorted):")
print(np.round(numerical_eigenvalues_sorted, 8))
print("\nAnalytically derived eigenvalues (sorted):")
print(np.round(analytical_eigenvalues, 8))

# Verify that the numerical and analytical results match
if np.allclose(analytical_eigenvalues, numerical_eigenvalues_sorted):
    print("\nVerification successful: Numerical and analytical eigenvalues match.\n")
else:
    print("\nVerification failed: Numerical and analytical eigenvalues do not match.\n")


print("--- 4. Rate of Relaxation ---")
print("The relaxation rate is determined by the second-largest eigenvalue.")

# Identify the eigenvalues
largest_eigenvalue = analytical_eigenvalues[0]
second_largest_eigenvalue = analytical_eigenvalues[1] # Note: lambda_1 = lambda_{N-1}

print(f"Largest eigenvalue (lambda_0): {largest_eigenvalue:.6f}")
print(f"Second-largest eigenvalue (lambda_1): {second_largest_eigenvalue:.6f}")
print("This value corresponds to n=1 in the formula lambda_n = cos(2*pi*n/N).\n")

print("The final equation for the second-largest eigenvalue is:")
print("lambda_1 = cos(2 * pi / N)")
print(f"For N = {N}, the calculation is:")
print(f"lambda_1 = cos(2 * {np.pi:.6f} / {N}) = cos({2 * np.pi / N:.6f}) = {np.cos(2 * np.pi / N):.6f}")
