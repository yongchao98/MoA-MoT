import numpy as np

def analyze_random_walk(N=20):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function:
    1. Constructs the transition probability matrix A.
    2. Verifies the analytic form of eigenvectors and eigenvalues.
    3. Calculates the rate of relaxation from the second-largest eigenvalue.
    """
    if N <= 2:
        print("N must be greater than 2 for a meaningful analysis.")
        return

    # --- Introduction ---
    print(f"--- Analysis for a Random Walk on a Circle with N={N} sites ---")
    print("The one-step transformation of the probability distribution P is P_t+1 = A @ P_t.")

    # --- Step 1: Construct the Transition Probability Matrix A ---
    print("\nStep 1: The Transition Probability Matrix A")
    # A[i, j] is the probability of moving from site j to site i.
    # In a symmetric random walk, from site j, one can move to (j-1)%N or (j+1)%N.
    # We use 0-indexing for sites: 0, 1, ..., N-1.
    A = np.zeros((N, N))
    for j in range(N):
        # Probability of moving right: j -> (j+1)%N
        A[(j + 1) % N, j] = 0.5
        # Probability of moving left: j -> (j-1)%N
        A[(j - 1) % N, j] = 0.5
    print("The transition probability matrix A has been constructed.")
    # print(A) # Uncomment to see the full matrix for small N

    # --- Step 2: Verify Eigenvectors and Eigenvalues ---
    print("\nStep 2: Verification of Eigenvectors and Eigenvalues")
    # The eigenvectors are v_n with components (v_n)_j = exp(i * k_n * j),
    # where k_n = 2 * pi * n / N, for n = 0, 1, ..., N-1.
    # The corresponding eigenvalue is lambda_n = cos(k_n).
    
    # We will verify this for the case n=1.
    n_test = 1
    k_n_test = 2 * np.pi * n_test / N
    lambda_n_analytic = np.cos(k_n_test)
    j_indices = np.arange(N)
    v_n_test = np.exp(1j * k_n_test * j_indices)

    # Check if A @ v_n = lambda_n * v_n
    lhs = A @ v_n_test
    rhs = lambda_n_analytic * v_n_test
    
    if np.allclose(lhs, rhs):
        print(f"Successfully verified that v_n with n={n_test} is an eigenvector.")
        print(f"The corresponding eigenvalue is cos(2*pi*{n_test}/{N}) = {lambda_n_analytic:.6f}")
    else:
        print(f"Verification failed for n={n_test}.")


    # --- Step 3: Calculate the Rate of Relaxation ---
    print("\nStep 3: Calculating the Rate of Relaxation")
    # The eigenvalues are lambda_n = cos(2*pi*n/N).
    # The largest eigenvalue (for n=0) is lambda_0 = cos(0) = 1.
    # The second-largest eigenvalue in magnitude corresponds to n=1 and n=N-1.
    n_sec = 1
    second_largest_eigenvalue = np.cos(2 * np.pi * n_sec / N)

    print(f"The largest eigenvalue is lambda_0 = cos(0) = 1.0")
    print("The second-largest eigenvalue (in magnitude) corresponds to n=1:")
    print(f"lambda_1 = cos(2 * pi * n / N)")
    print(f"lambda_1 = cos(2 * {np.pi:.6f} * {n_sec} / {N})")
    print(f"lambda_1 = {second_largest_eigenvalue:.6f}")

    # The rate of relaxation is determined by the spectral gap: 1 - |lambda_sec|.
    # Since lambda_1 is positive for N > 2, this is 1 - lambda_1.
    rate_of_relaxation = 1 - second_largest_eigenvalue

    print("\nThe rate of relaxation is given by 1 minus the second-largest eigenvalue:")
    print(f"Rate = 1 - lambda_1")
    print(f"Rate = 1 - {second_largest_eigenvalue:.6f}")
    print(f"Rate = {rate_of_relaxation:.6f}")
    
    # For comparison, the approximation for large N is 2*(pi/N)^2
    large_N_approx = 2 * (np.pi/N)**2
    print(f"\nFor large N, this rate can be approximated by 2*(pi/N)^2 = {large_N_approx:.6f}")
    
    return rate_of_relaxation

# --- Execution ---
# You can change the number of sites (N) here.
number_of_sites = 20
final_rate = analyze_random_walk(N=number_of_sites)

# The final answer in the required format.
# Note: The output below is the numerical result of the calculation above.
print(f"\n<<< {final_rate:.6f} >>>")
