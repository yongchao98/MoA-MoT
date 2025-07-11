import numpy as np

# Set the number of sites on the circle. You can change this value.
N = 5

# --- 1. One-Step Transformation ---
print("--- One-Step Transformation of Probability Distribution ---")
print(f"Consider a circle with N = {N} sites, indexed j = 0, 1, ..., {N-1}.")
print("At each step, the walker moves to an adjacent site with equal probability (0.5).")
print("Let p_j(t) be the probability of being at site j at time t.")
print("The one-step transformation equation is:")
# Output the equation with specific numbers
print(f"p_j(t+1) = 0.5 * p_((j-1) mod {N})(t) + 0.5 * p_((j+1) mod {N})(t)")
print("-" * 55)
print()

# --- 2. Transition Probability Matrix A ---
print("--- Transition Probability Matrix A ---")
# The matrix A is an N x N matrix where A[i, j] is the probability
# of transitioning from state j to state i.
A = np.zeros((N, N))
for j in range(N):
    # Probability of moving from j to j+1 (mod N)
    A[(j + 1) % N, j] = 0.5
    # Probability of moving from j to j-1 (mod N)
    A[(j - 1) % N, j] = 0.5

print(f"For N = {N}, the transition matrix A is:")
print(A)
print("-" * 55)
print()

# --- 3. Eigenvectors and Eigenvalues ---
print("--- Analytical Eigenvectors and Eigenvalues ---")
print("The eigenvectors v_n and eigenvalues lambda_n of this matrix can be found analytically.")
print("The proposed eigenvectors have components v_n(j) = exp(i * k_n * j),")
print(f"where i is the imaginary unit and k_n = 2 * pi * n / N for n = 0, ..., {N-1}.")
print("\nProof sketch:")
print("Applying the matrix A to an eigenvector v_n gives:")
print("  (A * v_n)_j = 0.5 * v_n(j-1) + 0.5 * v_n(j+1)")
print("              = 0.5 * [exp(i*k_n*(j-1)) + exp(i*k_n*(j+1))]")
print("              = 0.5 * exp(i*k_n*j) * [exp(-i*k_n) + exp(i*k_n)]")
print("Using Euler's formula, exp(ix) + exp(-ix) = 2*cos(x), this simplifies to:")
print("              = exp(i*k_n*j) * cos(k_n)")
print("              = cos(k_n) * v_n(j)")
print("\nThis shows that v_n is an eigenvector with eigenvalue lambda_n = cos(k_n).")
print(f"So, the eigenvalues are lambda_n = cos(2 * pi * n / {N}).")
print("-" * 55)
print()

# --- 4. Rate of Relaxation ---
print("--- Rate of Relaxation ---")
print("The convergence to the stationary distribution (uniform) is determined by the eigenvalues.")
print("The largest eigenvalue is lambda_0 = cos(0) = 1.")
print("The rate of relaxation is governed by the second-largest eigenvalue,")
print("as it corresponds to the slowest decaying mode.")
print(f"The second-largest eigenvalue occurs at n=1, so lambda_1 = cos(2*pi*1/{N}).")

# Calculate the second-largest eigenvalue
n_for_second_largest = 1
lambda_2nd_largest = np.cos(2 * np.pi * n_for_second_largest / N)

# The rate of relaxation is the "spectral gap", 1 - lambda_1.
relaxation_rate = 1 - lambda_2nd_largest

print("\nThe final equation for the relaxation rate (R) is R = 1 - lambda_1:")
# Output the equation with the calculated numbers
print(f"R = 1 - cos(2 * pi / {N}) = 1 - {lambda_2nd_largest:.5f} = {relaxation_rate:.5f}")
print("-" * 55)
