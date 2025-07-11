import numpy as np

# Set the number of sites on the circle.
# You can change this value to see how it affects the results.
N = 10

# --- Part 1: The One-Step Transformation and Transition Matrix A ---
print(f"--- Analysis of Random Walk on a Circle with N = {N} sites ---")
print("\nThe probability distribution p(t) transforms according to:")
print("p_j(t+1) = 0.5 * p_{j-1}(t) + 0.5 * p_{j+1}(t)  (indices are mod N)")
print("This can be written in matrix form as p(t+1) = A @ p(t),")
print("where A is the transition probability matrix.")

# Construct the transition matrix A
# A[i, j] represents the probability of moving from state j to state i.
# In our notation p(t+1) = A @ p(t), the rows of A are indexed by the new state.
A = np.zeros((N, N))
for i in range(N):
    # Probability of being at site i comes from site i-1 or i+1
    A[i, (i - 1) % N] = 0.5
    A[i, (i + 1) % N] = 0.5

print(f"\nThe transition matrix A for N={N} is:")
# Use suppress=True to make the output cleaner
with np.printoptions(precision=3, suppress=True):
    print(A)
print("-" * 60)

# --- Part 2: Verification of Eigenvectors and Eigenvalues ---
print("\n--- Verifying the Eigenvectors and Eigenvalues ---")
print("The analytical eigenvectors are v_n[j] = exp(i * 2*pi*n*j / N)")
print("The analytical eigenvalues are lambda_n = cos(2*pi*n / N)")

# We will verify this for n=1
n = 1
k_n = 2 * np.pi * n / N

# Analytical eigenvector v_n
v_n = np.exp(1j * k_n * np.arange(N))
# Analytical eigenvalue lambda_n
lambda_n = np.cos(k_n)

print(f"\nLet's test for n={n}:")
print(f"The equation to verify is: A @ v_{n} = lambda_{n} * v_{n}")
print(f"The analytical eigenvalue is lambda_{n} = cos(2*pi*{n}/{N}) = {lambda_n:.5f}")

# Calculate the left-hand side: A @ v_n
lhs = A @ v_n
# Calculate the right-hand side: lambda_n * v_n
rhs = lambda_n * v_n

print("\nLeft-hand side of the equation (A @ v_1):")
print("Each number is a component of the resulting vector (Real + Imaginary*j).")
for val in lhs:
    print(f"{val.real:9.5f} + {val.imag:9.5f}j")

print("\nRight-hand side of the equation (lambda_1 * v_1):")
print("Each number is a component of the resulting vector (Real + Imaginary*j).")
for val in rhs:
    print(f"{val.real:9.5f} + {val.imag:9.5f}j")

# Check if the results are numerically close
is_eigenvector = np.allclose(lhs, rhs)
print(f"\nAre the results equal (within tolerance)? {is_eigenvector}")
assert is_eigenvector, "Verification failed!"
print("The analytical form is correct.")
print("-" * 60)

# --- Part 3: Rate of Relaxation ---
print("\n--- Finding the Rate of Relaxation ---")
print("The largest eigenvalue is lambda_0 = cos(0) = 1.0, corresponding to the stationary state.")

# The rate of relaxation is determined by the second-largest eigenvalue.
second_largest_eigenvalue = np.cos(2 * np.pi * 1 / N)

print("The rate of relaxation is determined by the second-largest eigenvalue, which is lambda_1.")
print("\nFinal Equation for the rate-determining eigenvalue:")
print(f"lambda_1 = cos(2 * pi / N)")
print("\nCalculation for each number in the final equation:")
print(f"  pi ≈ {np.pi:.5f}")
print(f"  N  = {N}")
print(f"Result: lambda_1 = cos(2 * {np.pi:.5f} / {N}) = cos({2 * np.pi / N:.5f}) = {second_largest_eigenvalue:.5f}")
print("-" * 60)