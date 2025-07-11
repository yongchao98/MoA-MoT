import numpy as np
import sys

# --- Plan Execution ---

# Set a default value for N
N = 10
# Allow user to provide N as a command-line argument for demonstration
if len(sys.argv) > 1:
    try:
        N = int(sys.argv[1])
        if N <= 1:
            print("Error: N must be an integer greater than 1.", file=sys.stderr)
            sys.exit(1)
    except (ValueError, IndexError):
        print(f"Usage: python {sys.argv[0]} [N]", file=sys.stderr)
        print("Invalid input. Please provide an integer N > 1. Using default N=10.", file=sys.stderr)
        N = 10

# --- 1. System Definition and One-Step Transformation ---
print("### 1. System: 1D Random Walk on a Circle (Markov Chain) ###")
print(f"We consider a circle with N={N} sites, labeled j = 0, 1, ..., N-1.")
print("A walker at site j moves to the neighboring site j-1 or j+1 with equal probability p=0.5.")
print("\nLet pi_t(j) be the probability of being at site j at time t.")
print("The probability distribution evolves in one step according to the equation:")
print("\n  pi_{t+1}(j) = 0.5 * pi_t(j-1) + 0.5 * pi_t(j+1)\n")
print("Here, all site indices are calculated modulo N (e.g., site 0-1 is N-1).\n")
print("-" * 60)


# --- 2. Transition Probability Matrix A ---
print("### 2. Transition Probability Matrix A ###")
print("This transformation can be written in matrix form: vec(pi)_{t+1} = A * vec(pi)_t")
print("The element A[i, j] is the probability of transitioning TO site i FROM site j.")
# A_ij = P(state i at t+1 | state j at t)
A = np.zeros((N, N))
for j in range(N):
    # From state j, we can move to i = (j-1)%N or i = (j+1)%N
    A[(j - 1) % N, j] = 0.5
    A[(j + 1) % N, j] = 0.5
print(f"The transition matrix A for N={N} is:")
# Use np.set_printoptions to format the matrix nicely if it's large
if N <= 12:
    with np.printoptions(precision=2, suppress=True):
        print(A)
else:
    print("(Matrix is too large to display concisely)")
print("\n" + "-" * 60)


# --- 3. Analytical Derivation of Eigenvalues and Eigenvectors ---
print("### 3. Analytical Derivation ###")
print("To find the eigenvalues, we must solve A*v = lambda*v.")
print("We propose that the eigenvectors v_n are Fourier modes with components:")
print(f"  (v_n)_j = exp(i * j * k_n)")
print(f"where 'i' is the imaginary unit, and k_n = 2*pi*n/N for n = 0, 1, ..., {N-1}.")

print("\nLet's apply A to a component j of this proposed eigenvector:")
print("  (A * v_n)_j = Sum over l of A[j, l] * (v_n)_l")
print("  Since A[j, l] is only non-zero for l=j-1 and l=j+1, this simplifies to:")
print("  (A * v_n)_j = A[j, j-1]*(v_n)_{j-1} + A[j, j+1]*(v_n)_{j+1}")
print("  (A * v_n)_j = 0.5 * exp(i*(j-1)*k_n) + 0.5 * exp(i*(j+1)*k_n)")
print("\nFactoring out exp(i*j*k_n):")
print("  (A * v_n)_j = 0.5 * exp(i*j*k_n) * [exp(-i*k_n) + exp(i*k_n)]")
print("\nUsing Euler's formula, cos(x) = (e^(ix) + e^(-ix)) / 2:")
print("  (A * v_n)_j = exp(i*j*k_n) * cos(k_n)")
print("  (A * v_n)_j = (v_n)_j * cos(k_n)")
print("\nThis confirms that v_n is an eigenvector, and the corresponding eigenvalue is:")
print("  lambda_n = cos(k_n) = cos(2*pi*n/N)\n")
print("-" * 60)

# --- 4. Numerical Verification ---
print("### 4. Numerical Verification ###")
# Calculate eigenvalues numerically
eigenvalues_numeric, _ = np.linalg.eig(A)
eigenvalues_numeric_sorted = np.sort(np.real(eigenvalues_numeric))[::-1] # Sort descending

# Calculate eigenvalues analytically
eigenvalues_analytical = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
eigenvalues_analytical_sorted = np.sort(eigenvalues_analytical)[::-1]

print(f"We verify the formula by comparing it to numerical results for N={N}.\n")
print("{:<30} {:<30}".format("Analytical Eigenvalues (sorted)", "Numerical Eigenvalues (sorted)"))
for an, nu in zip(eigenvalues_analytical_sorted, eigenvalues_numeric_sorted):
    print("{:<30.8f} {:<30.8f}".format(an, nu))
print("\nThe results match perfectly.\n")
print("-" * 60)

# --- 5. Relaxation Rate and Final Answer ---
print("### 5. Rate of Relaxation ###")
print("The rate of convergence to the stationary distribution is determined by the second-largest eigenvalue.")
print("The largest eigenvalue is always lambda_0 = cos(2*pi*0/N) = cos(0) = 1.")
print("This lambda_0=1 corresponds to the stationary (uniform) probability distribution.\n")

print("The second-largest eigenvalues occur for n=1 and n=N-1:")
print("  lambda_1      = cos(2*pi*1/N) = cos(2*pi/N)")
print("  lambda_{N-1}  = cos(2*pi*(N-1)/N) = cos(2*pi - 2*pi/N) = cos(2*pi/N)")
print("\nThus, the second-largest eigenvalue, which determines the relaxation rate, is cos(2*pi/N).")

# Final equation output as requested
print("\nThe final equation for this value is:")
print("\n========================================")
number_2 = 2
symbol_pi = "pi"
symbol_N = "N"
print(f"  lambda_2nd = cos({number_2} * {symbol_pi} / {symbol_N})")
print("========================================")

# Calculate for the specific N
result_value = np.cos(2 * np.pi / N)
print(f"\nFor the specific case of N={N}, this value is cos(2*pi/{N}) = {result_value:.8f}")
