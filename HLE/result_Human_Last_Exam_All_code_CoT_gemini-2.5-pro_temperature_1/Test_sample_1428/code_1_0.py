import numpy as np

# Set the number of sites on the circle
N = 10
print(f"--- Random Walk on a Circle with N = {N} sites ---")

# --- Step 1: Construct the Transition Matrix A ---
print("\nStep 1: Constructing the Transition Matrix A")
# Initialize an N x N matrix with zeros
A = np.zeros((N, N))

# Populate the matrix based on the random walk rules
# A[i, j] is the probability of moving from state j to state i
for j in range(N):
    # From state j, can move to (j-1)%N or (j+1)%N with probability 0.5
    # So, we set the corresponding rows in column j to 0.5
    A[(j - 1) % N, j] = 0.5
    A[(j + 1) % N, j] = 0.5

print("The transition matrix A is:")
print(A)

# --- Step 2: Verify Eigenvectors and Eigenvalues Analytically ---
print("\nStep 2: Verifying the Eigenvector-Eigenvalue Relationship")
print("The eigenvectors v_n have components v_n,j = exp(i * 2*pi*n*j / N)")
print("The eigenvalues lambda_n are given by lambda_n = cos(2*pi*n / N)")

# Let's verify for n=1
n = 1
print(f"\nVerifying for n = {n}:")

# Proposed eigenvalue
k_n = (2 * np.pi * n) / N
lambda_n_analytical = np.cos(k_n)

# Proposed eigenvector
j_indices = np.arange(N)
v_n = np.exp(1j * k_n * j_indices)

# Calculate LHS: A * v_n
lhs = A @ v_n

# Calculate RHS: lambda_n * v_n
rhs = lambda_n_analytical * v_n

print(f"The analytical eigenvalue for n={n} is lambda_{n} = cos(2*pi*{n}/{N}) = {lambda_n_analytical:.6f}")
print("We will now check if A @ v_n equals lambda_n * v_n.")
# print("LHS (A @ v_n):", np.round(lhs, 6))
# print("RHS (lambda_n * v_n):", np.round(rhs, 6))
if np.allclose(lhs, rhs):
    print("Verification successful: A @ v_n is indeed equal to lambda_n * v_n.")
else:
    print("Verification failed.")

# --- Step 3: Find the Rate of Relaxation ---
print("\nStep 3: Finding the Rate of Relaxation")
print("The rate of relaxation is determined by the second-largest eigenvalue.")

# The eigenvalues are lambda_n = cos(2*pi*n/N) for n = 0, 1, ..., N-1
all_n = np.arange(N)
all_eigenvalues = np.cos(2 * np.pi * all_n / N)

# The largest eigenvalue corresponds to n=0
lambda_0 = all_eigenvalues[0]
print(f"The largest eigenvalue is for n=0: lambda_0 = cos(2*pi*0/{N}) = {lambda_0}")

# The second-largest eigenvalue corresponds to n=1 (and n=N-1)
# since cos(x) is a decreasing function for x in [0, pi].
n_second = 1
lambda_second_largest = np.cos(2 * np.pi * n_second / N)

print("\n--- Final Result ---")
print("The rate of relaxation is determined by the second-largest eigenvalue, which is for n=1:")
print(f"Rate = lambda_1 = cos(2 * pi * {n_second} / {N})")
print(f"Rate = cos({(2*n_second/N):.2f} * pi) = {lambda_second_largest}")

# For comparison, we can also find the eigenvalues numerically from the matrix A
numeric_eigenvalues, _ = np.linalg.eig(A)
numeric_eigenvalues.sort() # Sorts in ascending order
# The second largest is the second from the end of the sorted list
second_largest_numeric = numeric_eigenvalues[-2]
print(f"\nFor confirmation, the second-largest eigenvalue computed numerically from A is: {second_largest_numeric:.15f}")
print(f"The analytical value is: {lambda_second_largest:.15f}")
assert np.isclose(lambda_second_largest, second_largest_numeric)
print("The analytical result matches the numerical computation.")

<<<0.8090169943749475>>>