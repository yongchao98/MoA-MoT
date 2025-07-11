import numpy as np

# Set the number of sites for demonstration purposes.
# The following derivation is general for any integer N > 2.
N = 10

print("### Analysis of a 1D Random Walk on a Circle ###\n")

# Step 1: Define the model and the one-step transformation
print("--- 1. Model and Transformation ---")
print(f"We consider a one-dimensional random walk on a circle with N = {N} sites, indexed j = 0, 1, ..., N-1.")
print("At each step, the walker moves to an adjacent site (left or right) with equal probability of 1/2.")
print("Let P_t(j) be the probability of being at site j at time t.")
print("The one-step transformation of the probability distribution is:")
print("P_{t+1}(j) = (1/2) * P_t((j-1) mod N) + (1/2) * P_t((j+1) mod N)\n")

# Step 2: Define the transition matrix A
print("--- 2. Transition Probability Matrix A ---")
print("This process is a Markov chain described by the equation P_{t+1} = A * P_t, where P_t is the column vector of probabilities.")
print("The transition matrix A has elements A_ij = P(transition from site j to site i).")
print("A_ij = 1/2 if i and j are neighbors (i.e., i = (j +/- 1) mod N), and A_ij = 0 otherwise.")
# For demonstration, let's build and show the matrix for the chosen N
A = np.zeros((N, N))
for j in range(N):
    A[(j - 1) % N, j] = 0.5
    A[(j + 1) % N, j] = 0.5
print(f"For N = {N}, the matrix A is:\n{A}\n")

# Step 3: Find the eigenvectors and eigenvalues
print("--- 3. Eigenvectors and Eigenvalues of A ---")
print("We test if vectors v_n with components (v_n)_j = exp(i*j*k_n) are eigenvectors of A,")
print(f"where k_n = 2*pi*n/N for n = 0, 1, ..., N-1.")
print("\nLet's compute the j-th component of the product A * v_n:")
print("(A * v_n)_j = Sum over l [ A_jl * (v_n)_l ]")
print("          = (1/2) * (v_n)_(j-1) + (1/2) * (v_n)_(j+1)")
print("          = (1/2) * exp(i*(j-1)*k_n) + (1/2) * exp(i*(j+1)*k_n)")
print("          = exp(i*j*k_n) * [ (1/2) * (exp(-i*k_n) + exp(i*k_n)) ]")
print("\nUsing Euler's formula, cos(x) = (e^(ix) + e^(-ix)) / 2, the term in brackets is cos(k_n).")
print("(A * v_n)_j = exp(i*j*k_n) * cos(k_n) = (v_n)_j * cos(k_n)")
print("\nThis confirms that v_n is an eigenvector with the corresponding eigenvalue lambda_n = cos(k_n).")
print("The eigenvalues of the transition matrix A are therefore:")
print(f"lambda_n = cos(2 * pi * n / N) for n = 0, 1, ..., N-1\n")

# Step 4: Determine the rate of relaxation
print("--- 4. Rate of Relaxation ---")
print("The evolution of any initial probability distribution P_0 is given by P_t = A^t * P_0.")
print("In the eigenbasis, this becomes P_t = Sum over n [ c_n * (lambda_n)^t * v_n ].")
print("The largest eigenvalue is for n=0: lambda_0 = cos(0) = 1. This corresponds to the stationary distribution (uniform over all sites).")
print("All other eigenvalues have a magnitude |lambda_n| < 1 (for n > 0 and N > 2), so the corresponding terms decay to zero as t -> infinity.")
print("The rate of relaxation to the stationary distribution is governed by the term that decays the slowest.")
print("This corresponds to the eigenvalue with the largest magnitude less than 1, which is the second-largest eigenvalue.")
print("\nThe second-largest eigenvalue occurs for n=1 (and n=N-1):")
print("lambda_1 = cos(2 * pi * 1 / N)")
print("\nThis value determines the rate of relaxation. The closer lambda_1 is to 1, the slower the convergence.")
print("The final equation for the second-largest eigenvalue, which determines the relaxation rate, is:")

# Define numbers for the final equation as requested
number_two = 2
number_one = 1
pi_symbol = "pi"
print(f"lambda_1 = cos({number_two} * {pi_symbol} * {number_one} / N)")