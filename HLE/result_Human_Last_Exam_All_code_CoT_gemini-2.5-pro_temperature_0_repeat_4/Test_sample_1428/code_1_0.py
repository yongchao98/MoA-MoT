import numpy as np

def solve_random_walk_on_circle(N):
    """
    Analyzes the random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle. Must be greater than 1.
    """
    if not isinstance(N, int) or N <= 1:
        print("Error: N must be an integer greater than 1.")
        return

    # --- Step 1 & 2: Define the transformation and the transition matrix A ---
    # The probability distribution transforms as P_{t+1} = A * P_t.
    # For a site i, the probability is updated as P_{t+1}(i) = 0.5 * P_t(i-1) + 0.5 * P_t(i+1).
    # This means the matrix A has non-zero elements A[i, (i-1)%N] = 0.5 and A[i, (i+1)%N] = 0.5.
    # We use 0-based indexing for sites (0, 1, ..., N-1).
    
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1 + N) % N] = 0.5
        A[i, (i + 1) % N] = 0.5

    print("--- Task Analysis ---")
    print("1. The one-step transformation of the probability distribution P_t is P_{t+1} = A * P_t.")
    print(f"2. For N = {N}, the transition probability matrix A is:")
    print(A)
    print("\n")

    # --- Step 3: Verify the eigenvectors and find eigenvalues ---
    # The problem states the eigenvectors v_n have components (v_n)_j = exp(i*j*k_n)
    # with k_n = 2*pi*n/N. The corresponding eigenvalues are lambda_n = cos(k_n).
    # Let's verify this for a sample case, n=1.
    
    print("3. Verifying the eigenvectors and eigenvalues.")
    print("   The eigenvectors v_n have components (v_n)_j = exp(i*j*k_n) for k_n = 2*pi*n/N.")
    print("   The corresponding eigenvalues are lambda_n = cos(k_n).")
    
    n_test = 1
    k_n = 2 * np.pi * n_test / N
    
    # Construct the proposed eigenvector v_n
    v_n = np.array([np.exp(1j * j * k_n) for j in range(N)])
    
    # Calculate the left-hand side: A * v_n
    Av_n = A @ v_n
    
    # Calculate the right-hand side: lambda_n * v_n
    lambda_n = np.cos(k_n)
    lambda_v_n = lambda_n * v_n
    
    is_eigenvector = np.allclose(Av_n, lambda_v_n)
    
    print(f"\n   Verification for n = {n_test}:")
    print(f"   - Is v_{n_test} an eigenvector with eigenvalue lambda_{n_test}? {is_eigenvector}")
    if not is_eigenvector:
        print("     Verification failed. There might be a numerical precision issue or an error in the logic.")
    print("\n")

    # --- Step 4: Find the rate of relaxation ---
    # The rate of relaxation is determined by the second-largest eigenvalue.
    # The largest eigenvalue is lambda_0 = cos(0) = 1, corresponding to the stationary state.
    # The second-largest eigenvalue corresponds to n=1 (and n=N-1), which is lambda_1 = cos(2*pi/N).
    
    print("4. Finding the rate of relaxation.")
    print("   The rate is determined by the second-largest eigenvalue.")
    
    # Calculate all eigenvalues
    all_eigenvalues = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    
    # Find the unique eigenvalues and sort them in descending order
    unique_sorted_eigenvalues = np.sort(np.unique(all_eigenvalues))[::-1]
    
    lambda_max = unique_sorted_eigenvalues[0]
    second_lambda_max = unique_sorted_eigenvalues[1]

    print(f"\n   The largest eigenvalue is lambda_0 = cos(2*pi*0/{N}) = {lambda_max:.6f}")
    
    print("\n   The final equation for the second-largest eigenvalue is:")
    print(f"   lambda_1 = cos(2 * pi / N)")
    print(f"   For N = {N}, the value is:")
    print(f"   lambda_1 = cos(2 * {np.pi:.4f} / {N}) = {second_lambda_max:.6f}")

    print("\n--- Conclusion ---")
    print("The rate of relaxation towards the stationary distribution is determined by the")
    print(f"second-largest eigenvalue, which is {second_lambda_max:.6f} for N={N}.")
    print(f"The system's deviation from equilibrium decays proportionally to ({second_lambda_max:.6f})^t at each step t.")


if __name__ == '__main__':
    # You can change the number of sites (N) here
    number_of_sites = 10
    solve_random_walk_on_circle(number_of_sites)