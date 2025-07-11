import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes a 1D random walk on a circle, calculates its transition matrix,
    verifies its eigenvectors and eigenvalues, and finds the rate of relaxation.
    """
    # Set the number of sites on the circle
    N = 8

    # --- Introduction ---
    print(f"Analyzing a 1D random walk on a circle with N={N} sites.")
    print("-" * 60)

    # --- 1. One-step transformation ---
    print("1. One-Step Transformation of the Probability Distribution\n")
    print("Let p_i(t) be the probability of being at site i at time t.")
    print("From any site, the walker moves to an adjacent site with equal probability.")
    print("The one-step transformation equation is:")
    # Output the numbers in the final equation
    print("p_i(t+1) = (1/2) * p_{i-1}(t) + (1/2) * p_{i+1}(t)\n")
    print("-" * 60)

    # --- 2. Transition Probability Matrix A ---
    print("2. Transition Probability Matrix A\n")
    print("This transformation can be written in matrix form: p(t+1) = A * p(t).")
    print(f"The {N}x{N} transition matrix A for N={N} is constructed as follows:\n")
    # Create the transition matrix A
    A = np.zeros((N, N))
    for i in range(N):
        # A_ij is the probability of moving from j to i.
        # So, from column j, we can go to (j-1) or (j+1).
        A[(i - 1 + N) % N, i] = 0.5 # Moving from i to i-1
        A[(i + 1) % N, i] = 0.5       # Moving from i to i+1
    print(A)
    print("\n" + "-" * 60)

    # --- 3. Eigenvectors and Eigenvalues ---
    print("3. Verifying Eigenvectors and Eigenvalues\n")
    print("The eigenvectors v_n have components v_n_j = exp(i * k_n * j)")
    print("with k_n = 2*pi*n/N for n = 0, ..., N-1.")
    print("The corresponding analytical eigenvalue is lambda_n = cos(2*pi*n/N).\n")
    print("Verifying the eigenvalue equation A * v_n = lambda_n * v_n:")

    j_indices = np.arange(N) # site indices
    for n in range(N):
        k_n = 2 * np.pi * n / N
        # Theoretical eigenvector and eigenvalue
        v_n = np.exp(1j * k_n * j_indices)
        lambda_n_analytical = np.cos(k_n)
        
        # Check if A * v_n = lambda_n * v_n
        lhs = A @ v_n
        rhs = lambda_n_analytical * v_n
        
        is_eigen_pair = np.allclose(lhs, rhs)
        print(f"\nFor n={n}:")
        print(f"  Eigenvalue equation: lambda_{n} = cos(2*pi*{n}/{N})")
        print(f"  Calculated analytical eigenvalue: {lambda_n_analytical:.6f}")
        print(f"  Verification check (A*v == lambda*v): {'Passed' if is_eigen_pair else 'Failed'}")
        
    print("\n" + "-" * 60)

    # --- 4. Rate of Relaxation ---
    print("4. Rate of Relaxation\n")
    print("The rate of relaxation towards the stationary distribution is")
    print("determined by the second-largest eigenvalue of the matrix A.")
    
    # Calculate all analytical eigenvalues
    eigenvalues = np.cos(2 * np.pi * np.arange(N) / N)
    
    # The largest eigenvalue (lambda_0 for n=0) corresponds to the stationary state
    lambda_max = eigenvalues[0]
    print(f"\nThe largest eigenvalue (for n=0) is: lambda_0 = cos(2*pi*0/{N}) = {lambda_max:.6f}")

    # The second-largest eigenvalue is for n=1 (and n=N-1)
    unique_sorted_eigenvalues = np.sort(np.unique(eigenvalues))[::-1]
    lambda_second_largest = unique_sorted_eigenvalues[1]

    print("The second-largest eigenvalue corresponds to the mode n=1.")
    print("The equation for this eigenvalue is: lambda_1 = cos(2*pi*1/N)")
    print(f"For N={N}, its value is: cos(2*pi*1/{N}) = {lambda_second_largest:.6f}")

if __name__ == '__main__':
    solve_random_walk_on_circle()

<<<0.707107>>>