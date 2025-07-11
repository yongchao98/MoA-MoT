import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes a 1D random walk on a circle, calculates its transition matrix,
    verifies its eigenvalues, and finds the rate of relaxation.
    """
    # Define the number of sites on the circle
    N = 10
    
    print(f"Analyzing a 1D random walk on a circle with N = {N} sites.")
    print("-" * 60)

    # --- Step 1: One-step transformation ---
    print("Step 1: The One-Step Transformation")
    print("Let p(t) be the vector of probabilities of being at each site at time t.")
    print("The system evolves according to p(t+1) = A @ p(t), where A is the")
    print("transition probability matrix.")
    print("\nFor any site i, the probability of being there at t+1 is the sum of")
    print("probabilities of coming from its neighbors (i-1 and i+1) at time t.")
    print("p_i(t+1) = 0.5 * p_{i-1}(t) + 0.5 * p_{i+1}(t)")
    print("-" * 60)

    # --- Step 2: Construct the transition matrix A ---
    # We use 0-based indexing for sites (0, 1, ..., N-1) for convenience.
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5
        
    print("Step 2: The Transition Probability Matrix A")
    print("The matrix A where A_ij is the probability of moving from j to i in one step:")
    print(A)
    print("-" * 60)

    # --- Step 3: Eigenvectors and Eigenvalues ---
    print("Step 3: Theoretical Eigenvectors and Eigenvalues")
    print("The eigenvectors of A are v_n, with components (v_n)_j = exp(i * k_n * j)")
    print("where k_n = 2 * pi * n / N for n = 0, 1, ..., N-1.")
    print("\nThe corresponding eigenvalues are lambda_n = cos(k_n) = cos(2 * pi * n / N).")
    print("-" * 60)

    # --- Step 4: Verification for n=1 ---
    print("Step 4: Verification for a non-trivial eigenvector (n=1)")
    n = 1
    # Wave number k_n for n=1
    k_n = 2 * np.pi * n / N
    
    # Theoretical eigenvector v_n
    v_n = np.exp(1j * k_n * np.arange(N))
    
    # Theoretical eigenvalue lambda_n
    lambda_n_theoretical = np.cos(k_n)
    
    # Action of A on v_n
    action_of_A = A @ v_n
    
    # Expected result
    expected_result = lambda_n_theoretical * v_n
    
    # Check if they are equal
    is_eigenvector = np.allclose(action_of_A, expected_result)
    
    print(f"For n={n}, we check if A @ v_{n} = lambda_{n} * v_{n}")
    print(f"Theoretical eigenvalue lambda_{n} = {lambda_n_theoretical:.6f}")
    print(f"Is v_{n} an eigenvector? {is_eigenvector}")
    print("-" * 60)

    # --- Step 5: Rate of Relaxation ---
    print("Step 5: Finding the Rate of Relaxation")
    print("The relaxation to the stationary state (uniform distribution) is governed")
    print("by the second-largest eigenvalue, as it represents the slowest-decaying mode.")
    
    # Largest eigenvalue (n=0)
    lambda_max = np.cos(2 * np.pi * 0 / N)
    print(f"\nThe largest eigenvalue is lambda_0 = cos(0) = {lambda_max}")

    # Second-largest eigenvalue (n=1 or n=N-1)
    lambda_second_largest = np.cos(2 * np.pi * 1 / N)
    
    print("The second-largest eigenvalue corresponds to n=1 (and n=N-1).")
    print("\nFinal equation and calculation for the relaxation rate determinant:")
    print(f"lambda_1 = cos(2 * pi / N)")
    print(f"         = cos(2 * {np.pi:.6f} / {N})")
    print(f"         = cos({2 * np.pi / N:.6f})")
    print(f"         = {lambda_second_largest:.6f}")
    print("\nThis second-largest eigenvalue determines the rate of relaxation.")
    print("-" * 60)

if __name__ == "__main__":
    solve_random_walk_on_circle()
