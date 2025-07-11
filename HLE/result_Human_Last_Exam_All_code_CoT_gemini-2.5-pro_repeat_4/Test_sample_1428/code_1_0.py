import numpy as np

def random_walk_analysis(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    print(f"Analyzing a 1D random walk on a circle with N = {N} sites.")
    print("-" * 50)

    # 1. One-step transformation
    print("1. One-Step Transformation of the Probability Distribution")
    print("The probability P_j at site j evolves according to:")
    print("  P_j(t+1) = 0.5 * P_{(j-1) mod N}(t) + 0.5 * P_{(j+1) mod N}(t)")
    print("In matrix form, this is p(t+1) = A @ p(t), where p is the probability vector.\n")

    # 2. Transition Probability Matrix A
    print("2. Transition Probability Matrix A")
    # A_ij is the probability of moving from state j to state i.
    A = np.zeros((N, N))
    for j in range(N):
        # From site j, can move to (j+1)%N or (j-1)%N
        # So, column j will have 0.5 at rows (j+1)%N and (j-1)%N
        A[(j + 1) % N, j] = 0.5
        A[(j - 1) % N, j] = 0.5
    print("The transition matrix A is:")
    print(A)
    print("\n")

    # 3. Eigenvectors and Eigenvalues
    print("3. Eigenvectors and Eigenvalues")
    print("The eigenvectors v_n of this matrix have components (v_n)_j given by:")
    print(f"  (v_n)_j = exp(i * k_n * j)   for j in [0, ..., {N-1}]")
    print(f"where k_n = 2 * pi * n / N   for n in [0, ..., {N-1}].")
    print("\nThe action of A on v_n shows that the corresponding eigenvalue is:")
    print("  lambda_n = cos(k_n) = cos(2 * pi * n / N)")
    print("\nLet's verify this numerically for n=1:")
    n = 1
    k_n = 2 * np.pi * n / N
    lambda_n = np.cos(k_n)
    v_n = np.array([np.exp(1j * k_n * j) for j in range(N)])
    
    lhs = A @ v_n
    rhs = lambda_n * v_n
    
    print(f"  For n={n}, the eigenvalue lambda_1 is cos(2*pi/{N}) = {lambda_n:.6f}")
    print(f"  Is A @ v_1 == lambda_1 * v_1?  ->  {np.allclose(lhs, rhs)}")
    print("\n")

    # 4. Rate of Relaxation
    print("4. Rate of Relaxation")
    print("The rate of relaxation to the stationary distribution is determined by the")
    print("second-largest eigenvalue.")
    print("The largest eigenvalue is lambda_0 = cos(0) = 1.")
    print("The second-largest eigenvalue corresponds to n=1 (and n=N-1), which gives")
    print("the slowest-decaying non-stationary mode.\n")

    # Final Equation Output
    print("Final Equation for the second-largest eigenvalue:")
    two = 2
    pi_val = np.pi
    print(f"  lambda_1 = cos({two} * {pi_val:.6f} / {N})")

    second_largest_eigenvalue_val = np.cos(2 * np.pi / N)
    print(f"\nFor N={N}, the value of this eigenvalue is: {second_largest_eigenvalue_val:.6f}")


# Execute the analysis for a chosen N
random_walk_analysis(N=10)