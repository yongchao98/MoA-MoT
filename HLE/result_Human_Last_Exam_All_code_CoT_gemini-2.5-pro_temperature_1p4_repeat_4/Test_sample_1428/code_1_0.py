import numpy as np

# --- Configuration ---
# Number of sites on the circle.
# We will use a concrete value to demonstrate the concepts numerically.
N = 10

# --- Introduction ---
print("This script analyzes a 1D symmetric random walk on a circle.")
print(f"We consider a circle with N = {N} sites, indexed j = 0, 1, ..., {N-1}.")
print("At each time step, a walker at site j moves to site (j-1)%N or (j+1)%N, each with probability 1/2.")
print("-" * 60)

# --- Step 1: One-step transformation of probability distribution ---
print("Step 1: The one-step transformation of the probability distribution\n")
print("Let P_j(t) be the probability of the walker being at site j at time t.")
print("The probability at the next time step, t+1, is given by:")
print(f"  P_j(t+1) = 0.5 * P_[(j-1)%{N}](t) + 0.5 * P_[(j+1)%{N}](t)")
print("\nThis equation describes how the probability distribution evolves in a single step.")
print("-" * 60)

# --- Step 2: The transition probability matrix A ---
print("Step 2: The transition probability matrix A\n")
print("The transformation can be expressed as a matrix equation: P(t+1) = A * P(t),")
print("where P(t) is the column vector of probabilities [P_0(t), P_1(t), ..., P_{N-1}(t)].")
print("The matrix element A[i, j] is the probability of transitioning FROM site j TO site i.")
print("So, A[i, j] is 0.5 if i = (j-1)%N or i = (j+1)%N, and 0 otherwise.")

# Construct the transition matrix A
# This corresponds to the equation P_i(t+1) = sum_j A_ij P_j(t).
A = np.zeros((N, N))
for i in range(N):
    # A walker at site (i-1)%N can move to site i.
    A[i, (i - 1) % N] = 0.5
    # A walker at site (i+1)%N can move to site i.
    A[i, (i + 1) % N] = 0.5

print(f"\nFor N = {N}, the transition matrix A is:")
# Setting print options for better display
np.set_printoptions(precision=2, suppress=True)
print(A)
print("-" * 60)

# --- Step 3: Verify the eigenvectors and eigenvalues ---
print("Step 3: Verifying the eigenvectors and eigenvalues\n")
print("It can be shown analytically that the eigenvectors and eigenvalues of A are:")
print(f"Eigenvector v_n corresponding to site j: (v_n)_j = exp(i * j * 2*pi*n/{N})")
print(f"Eigenvalue lambda_n: lambda_n = cos(2*pi*n/{N})")
print(f"for n = 0, 1, ..., {N-1}.\n")

print("Let's numerically verify this for all n:")
verification_passed = True
for n in range(N):
    # Theoretical eigenvector and eigenvalue
    k_n = 2 * np.pi * n / N
    lambda_n_theory = np.cos(k_n)
    v_n_theory = np.exp(1j * np.arange(N) * k_n)
    
    # Apply the matrix A to the theoretical eigenvector
    A_v = A @ v_n_theory
    
    # Compare with lambda * v
    lambda_v = lambda_n_theory * v_n_theory
    
    # Check if they are close (due to floating point arithmetic)
    if not np.allclose(A_v, lambda_v):
        print(f"Verification FAILED for n={n}")
        verification_passed = False
        break

if verification_passed:
    print("Verification successful: A @ v_n = lambda_n * v_n for all n.")
    print("This confirms the analytical solution.")
print("-" * 60)

# --- Step 4: Find the rate of relaxation ---
print("Step 4: Finding the rate of relaxation\n")
print("The convergence to the stationary distribution (uniform) is governed by the second-largest eigenvalue in magnitude.")
print("The stationary distribution corresponds to the largest eigenvalue, lambda_0 = cos(0) = 1.")

# The second-largest eigenvalue corresponds to n=1 (and n=N-1)
# which are equal: cos(2*pi/N) == cos(2*pi*(N-1)/N)
second_largest_eigenvalue = np.cos(2 * np.pi * 1 / N)
n_sl = 1

print("\nThe second-largest eigenvalue is for n=1 (and n=N-1):")
print(f"  lambda_{n_sl} = cos(2 * pi * {n_sl} / {N})")
print(f"  lambda_{n_sl} = cos({(2 * np.pi / N):.4f})")
print(f"  lambda_{n_sl} = {second_largest_eigenvalue:.4f}")

print("\nThe rate of relaxation is defined by the spectral gap: 1 - |lambda_SL|.")
relaxation_rate = 1 - second_largest_eigenvalue

print("\nThe relaxation rate is therefore:")
print(f"  Rate = 1 - lambda_{n_sl}")
print(f"  Rate = 1 - {second_largest_eigenvalue:.4f}")
print(f"  Rate = {relaxation_rate:.4f}")
print("-" * 60)

# Final answer formatted as requested
final_answer = 1 - np.cos(2 * np.pi / N)
print(f"<<<{final_answer}>>>")