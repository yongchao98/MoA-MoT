import numpy as np

# You can change the number of sites N here
N = 50

# --- Part 1: Define and Print the Transition Matrix A ---
# A_ij = P(next=i | current=j)
# So A_ij is 1/2 if i = (j-1)%N or i = (j+1)%N, and 0 otherwise.
A = np.zeros((N, N))
for j in range(N):
    # Probability 1/2 to move to the site on the "left"
    A[(j - 1) % N, j] = 0.5
    # Probability 1/2 to move to the site on the "right"
    A[(j + 1) % N, j] = 0.5

print(f"--- Transition Matrix A for N={N} (showing top-left 10x10) ---")
print(A[:10, :10])
print("...\n")

# --- Part 2: Verify the Eigenvectors and Eigenvalues ---
# We show that v_n with components (v_n)_j = exp(i*j*k_n) is an eigenvector
# of A with eigenvalue lambda_n = cos(k_n), where k_n = 2*pi*n/N.
# We will verify this numerically for the case n=1.

n_test = 1
k_n = 2 * np.pi * n_test / N

print(f"--- Verifying Eigenvector/Eigenvalue for n={n_test} ---")

# The proposed eigenvalue from the formula
lambda_n_analytic = np.cos(k_n)

# The proposed eigenvector from the formula
j_indices = np.arange(N)
v_n = np.exp(1j * j_indices * k_n)

# Apply the matrix A to the eigenvector: A * v_n
A_v = A @ v_n

# Multiply the eigenvector by the eigenvalue: lambda_n * v_n
lambda_v = lambda_n_analytic * v_n

# Check if the two results are the same
is_eigenvector = np.allclose(A_v, lambda_v)

print(f"Proposed eigenvalue lambda_{n_test} = cos(2*pi*{n_test}/{N}) = {lambda_n_analytic:.6f}")
print(f"We check if A @ v_n is equal to lambda_n * v_n.")
print(f"Are the results numerically equal? {is_eigenvector}\n")


# --- Part 3: Find the Rate of Relaxation ---
# The rate of relaxation is determined by the second-largest eigenvalue.
# This corresponds to n=1 (and n=N-1), which is lambda_1 = cos(2*pi/N).

# Final equation for the rate of relaxation
rate_of_relaxation = np.cos(2 * np.pi / N)

print("--- Rate of Relaxation (Second-Largest Eigenvalue) ---")
print("The rate is given by the equation: lambda_sl = cos(2 * pi / N)")
print(f"For N = {N}:")
print(f"lambda_sl = cos(2 * pi / {N})")
print(f"The calculated value is: {rate_of_relaxation}")

# Final numeric answer as requested by the prompt format
#<<<"%.10f" % rate_of_relaxation>>>