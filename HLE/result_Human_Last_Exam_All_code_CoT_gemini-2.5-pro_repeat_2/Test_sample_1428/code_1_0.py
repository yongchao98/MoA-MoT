import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes a 1D random walk on a circle with N sites.
    This script demonstrates the solution for a fixed N=5.
    """
    # Define the number of sites on the circle
    N = 5

    print(f"Solving the random walk on a circle problem for N = {N} sites.")
    print("=" * 60)

    # Step 1: Write the one-step transformation of the probability distribution
    print("Step 1: The one-step transformation of the probability distribution")
    print("Let P_t(i) be the probability of being at site i at time t.")
    print("For a symmetric random walk (probability 0.5 to move left, 0.5 to move right),")
    print("the probability at time t+1 is given by:")
    print(f"  P_{{t+1}}(i) = 0.5 * P_t(i-1) + 0.5 * P_t(i+1)")
    print("where the site indices are taken modulo N.")
    print("\nIn vector form, this is P_{t+1} = A * P_t, where A is the transition matrix.")
    print("-" * 60)

    # Step 2: Compute the transition probability matrix A
    print("Step 2: The transition probability matrix A")
    print(f"For N={N}, the matrix A is a {N}x{N} matrix where A_ij = 0.5 if site i")
    print("is a neighbor of site j, and 0 otherwise.")
    
    # We use 0-based indexing for arrays: states are 0, 1, ..., N-1
    A = np.zeros((N, N))
    for j in range(N):  # For each 'from' state j (matrix column)
        # The 'to' states i are the previous and next sites
        i_prev = (j - 1 + N) % N
        i_next = (j + 1) % N
        # The probability of moving from j to i_prev is 0.5
        A[i_prev, j] = 0.5
        # The probability of moving from j to i_next is 0.5
        A[i_next, j] = 0.5
        
    print("The calculated matrix A is:")
    print(A)
    print("-" * 60)

    # Step 3: Show that the eigenvectors are e^{ijk_n} and find eigenvalues
    print("Step 3: Eigenvectors and Eigenvalues of A")
    print("The transition matrix A is a circulant matrix. Its eigenvectors v_n are known")
    print("to be the vectors of the Discrete Fourier Transform.")
    print(f"Assuming 'l' in the prompt is the imaginary unit 'i', the j-th component")
    print(f"of the n-th eigenvector (for j,n = 0..{N-1}) is:")
    print(f"  (v_n)_j = exp(i * j * k_n), where k_n = 2*pi*n/N.")
    print("\nTo confirm this and find the eigenvalues, we apply A to v_n.")
    print("The j-th component of the vector (A * v_n) is:")
    print("  (A*v_n)_j = sum_{l=0..N-1} A_jl * (v_n)_l")
    print("            = 0.5 * (v_n)_{j-1} + 0.5 * (v_n)_{j+1}  (indices mod N)")
    print("            = 0.5 * exp(i * (j-1) * k_n) + 0.5 * exp(i * (j+1) * k_n)")
    print("            = exp(i * j * k_n) * 0.5 * [exp(-i * k_n) + exp(i * k_n)]")
    print("Using Euler's formula, cos(x) = (e^{ix} + e^{-ix})/2, this simplifies to:")
    print("  (A*v_n)_j = exp(i * j * k_n) * cos(k_n)")
    print("            = (v_n)_j * cos(2*pi*n/N)")
    print("\nThis confirms that v_n is an eigenvector with the corresponding eigenvalue lambda_n:")
    print("  lambda_n = cos(2 * pi * n / N)")
    print("-" * 60)

    # Step 4: Find the rate of relaxation
    print("Step 4: The rate of relaxation")
    print("The convergence to the stationary distribution (which corresponds to the largest")
    print("eigenvalue, lambda_0 = 1) is governed by the second-largest eigenvalue.")
    print("The largest eigenvalue is lambda_0 = cos(0) = 1.")
    print("The second-largest eigenvalue corresponds to n=1 (and n=N-1), which is cos(2*pi/N).")
    print("This eigenvalue determines the rate of relaxation.")
    
    print("\nThe final equation for the second-largest eigenvalue is:")
    n = 1
    pi_val = np.pi
    print(f"  lambda_1 = cos(2 * pi * n / N)")
    
    print(f"\nPlugging in the numbers for our example:")
    print(f"  n = {n}")
    print(f"  pi = {pi_val}")
    print(f"  N = {N}")

    second_largest_eigenvalue = np.cos(2 * pi_val * n / N)

    print(f"\n  lambda_1 = cos(2 * {pi_val:.5f} * {n} / {N}) = cos({2 * pi_val / N:.5f}) = {second_largest_eigenvalue:.5f}")
    print("\nThis value represents the relaxation rate of the slowest decaying mode.")

# Execute the function
solve_random_walk_on_circle()