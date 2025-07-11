import numpy as np

# Set the number of sites on the circle. You can change this value.
N = 10

# --- Introduction ---
print("--- 1D Random Walk on a Circle Analysis ---")
print(f"This script analyzes the system for N = {N} sites.")
print("The one-step probability transformation is: p_j(t+1) = 0.5 * p_{j-1}(t) + 0.5 * p_{j+1}(t)\n")

# --- 1. Compute and display the transition matrix A ---
print("1. The transition probability matrix A:")
# Initialize an N x N zero matrix
A = np.zeros((N, N))
# Populate the matrix based on the random walk rule
for j in range(N):
    # Transition from j-1 to j
    A[j, (j - 1 + N) % N] = 0.5
    # Transition from j+1 to j
    A[j, (j + 1 + N) % N] = 0.5
print(A)

# --- 2. Verify the eigenvectors and eigenvalues ---
print("\n2. Verifying the eigenvectors v_n and eigenvalues lambda_n = cos(2*pi*n/N):")
print("   We check if A @ v_n is numerically close to lambda_n * v_n for each n.")
all_verified = True
for n in range(N):
    k_n = 2 * np.pi * n / N
    # The predicted eigenvalue from the formula
    lambda_n_predicted = np.cos(k_n)

    # The predicted eigenvector from the formula
    v_n = np.array([np.exp(1j * j * k_n) for j in range(N)])

    # Numerically check if A @ v_n = lambda_n * v_n
    is_eigenvector = np.allclose(A @ v_n, lambda_n_predicted * v_n)
    if not is_eigenvector:
        all_verified = False
        
    print(f"   For n={n}: lambda_{n} = {lambda_n_predicted:.4f}. Verified: {is_eigenvector}")

if all_verified:
    print("\nVerification successful: The derived formula for eigenvectors and eigenvalues holds true.")
else:
    print("\nVerification failed for one or more eigenvectors.")


# --- 3. Find the rate of relaxation ---
print("\n3. Calculating the rate of relaxation:")
print("The rate is determined by the second-largest eigenvalue, λ₁ = cos(2*π/N).")

# Values used in the final calculation
two = 2
pi_val = np.pi

# Displaying the final equation with its numbers
print("\nFinal Equation Breakdown:")
print(f"For N = {N}:")
print(f"λ₁ = cos({two} * π / {N})")
print(f"λ₁ = cos({two} * {pi_val:.4f} / {N})")
print(f"λ₁ = cos({two * pi_val / N:.4f})")
final_value = np.cos(two * pi_val / N)
print(f"λ₁ = {final_value:.4f}")

<<<cos(2*pi/N)>>>