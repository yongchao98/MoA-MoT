import numpy as np

# Set the number of sites on the circle. You can change this value.
N = 10

# --- Part 1: One-Step Transformation Equation ---
print("The one-step transformation of the probability distribution pi over the sites is:")
print("pi_j(t+1) = 0.5 * pi_{j-1}(t) + 0.5 * pi_{j+1}(t), with indices modulo N.")
print("-" * 40)

# --- Part 2: Compute the Transition Probability Matrix A ---
# A[i, j] = Probability of transitioning from site j to site i.
A = np.zeros((N, N))
for j in range(N):
    # Walker at j moves to j-1 (i=j-1) or j+1 (i=j+1) with probability 0.5
    A[(j - 1) % N, j] = 0.5
    A[(j + 1) % N, j] = 0.5

print(f"Transition Matrix A for N = {N}:")
print(A)
print("-" * 40)

# --- Part 3: Verify Eigenvectors and Eigenvalues ---
# Let's numerically verify the relationship A*v_n = lambda_n*v_n for n=1.
# The claimed eigenvector v_n has components (v_n)_j = exp(i * j * k_n)
# The claimed eigenvalue is lambda_n = cos(k_n), where k_n = 2*pi*n/N.

n_to_verify = 1
k_n = 2 * np.pi * n_to_verify / N
j_indices = np.arange(N) # Represents sites 0, 1, ..., N-1

# Construct the theoretical eigenvector v_n
v_n = np.exp(1j * j_indices * k_n)

# Calculate the theoretical eigenvalue lambda_n
lambda_n = np.cos(k_n)

# Calculate the left side of the equation: A @ v_n
left_side = A @ v_n

# Calculate the right side of the equation: lambda_n * v_n
right_side = lambda_n * v_n

print(f"Verifying the eigenvector equation for n = {n_to_verify}: A*v_n = lambda_n*v_n")
# We use np.allclose for comparing floating-point complex numbers
is_eigenvector = np.allclose(left_side, right_side)
print(f"Verification successful: {is_eigenvector}")
# The check `np.allclose` confirms that the two vectors are equal within a small tolerance.
print("-" * 40)


# --- Part 4: Find the Rate of Relaxation ---
# The eigenvalues are lambda_n = cos(2*pi*n/N) for n = 0, ..., N-1.
eigenvalues_theoretical = np.cos(2 * np.pi * np.arange(N) / N)

# The largest eigenvalue corresponds to n=0: lambda_0 = cos(0) = 1.
largest_eigenvalue = eigenvalues_theoretical[0]

# The rate of relaxation is determined by the second-largest eigenvalue, which corresponds to n=1.
second_largest_eigenvalue = np.cos(2 * np.pi / N)

print("The rate of relaxation is determined by the second-largest eigenvalue.")
print(f"Largest eigenvalue (for n=0): {largest_eigenvalue:.4f}")
print("Final equation for the second-largest eigenvalue: lambda_1 = cos(2 * pi / N)")
print(f"For N = {N}, the equation is: lambda_1 = cos(2 * {np.pi:.4f} / {N})")
print(f"The resulting value is: {second_largest_eigenvalue:.4f}")
