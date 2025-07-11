import numpy as np

def analyze_random_walk_on_circle(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function constructs the transition matrix, verifies its eigenvectors
    and eigenvalues, and calculates the rate of relaxation determined by the
    second-largest eigenvalue.

    Args:
        N (int): The number of sites on the circle.
    """
    if not isinstance(N, int) or N < 3:
        print("Error: N must be an integer of 3 or greater.")
        return

    # Set print options for better readability of complex numbers and floats
    np.set_printoptions(precision=5, suppress=True)

    print(f"--- Analysis for a 1D Random Walk on a Circle with N={N} sites ---")
    print("\nStep 1: One-Step Transformation and Transition Matrix A")
    print("The probability distribution P evolves as P_{t+1} = A * P_t.")
    print("For a symmetric walk on a circle, a particle at site j moves to (j-1)%N or (j+1)%N with probability 0.5 each.")
    print("The element A[i, j] is the probability of moving from site j to site i.")

    # Construct the transition matrix A
    # A[i, j] = P(to i | from j)
    A = np.zeros((N, N))
    for j in range(N):
        A[(j - 1) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5

    print("\nThe transition matrix A is:")
    # Display smaller matrices, for larger ones just show the shape
    if N <= 10:
        print(A)
    else:
        print(f"<An {N}x{N} matrix, too large to display>")


    print("\nStep 2: Verify Eigenvectors and Eigenvalues")
    print("The eigenvectors v_n have components v_n[j] = exp(i*j*k_n) where k_n = 2*pi*n/N.")
    print("The corresponding eigenvalues are lambda_n = cos(k_n).")
    print("Let's verify this for n=1:")

    n = 1
    k_n = 2 * np.pi * n / N
    lambda_n_theory = np.cos(k_n)

    # Construct the theoretical eigenvector v_n
    j_indices = np.arange(N)
    v_n = np.exp(1j * j_indices * k_n)

    # Perform the matrix-vector multiplication A * v_n
    Av = A @ v_n
    # Calculate the theoretical result lambda_n * v_n
    lambda_v = lambda_n_theory * v_n

    print(f"\nFor n={n}, k_1 = 2*pi/{N} = {k_n:.5f}")
    print(f"The theoretical eigenvalue is lambda_1 = cos({k_n:.5f}) = {lambda_n_theory:.5f}")
    print("\nWe will now check if the equation A * v_1 = lambda_1 * v_1 holds.")
    
    # To avoid printing long vectors, we can check the result programmatically
    if np.allclose(Av, lambda_v):
        print("Verification successful: A * v_1 is indeed equal to lambda_1 * v_1.")
    else:
        print("Verification failed.")
        print("A * v_1 = ", Av)
        print("lambda_1 * v_1 = ", lambda_v)

    print("\nStep 3: Find the Rate of Relaxation")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    print("The eigenvalues are lambda_n = cos(2*pi*n/N) for n = 0, 1, ..., N-1.")
    
    # The largest eigenvalue is for n=0: lambda_0 = cos(0) = 1.
    # The second-largest eigenvalue is for n=1 (and n=N-1): lambda_1 = cos(2*pi/N).
    second_largest_eigenvalue = np.cos(2 * np.pi / N)

    print("\nThe largest eigenvalue is lambda_0 = cos(0) = 1.0, which corresponds to the stationary state.")
    print("The second-largest eigenvalue determines the relaxation rate.")

    print("\nFinal Equation and Result:")
    print(f"Second-largest eigenvalue = cos(2 * pi / N)")
    print(f"                        = cos(2 * {np.pi:.5f} / {N})")
    print(f"                        = cos({2 * np.pi / N:.5f})")
    print(f"                        = {second_largest_eigenvalue:.5f}")
    print("\nThis value determines the rate of relaxation of the system to its stationary state.")


if __name__ == '__main__':
    # You can change the number of sites (N) here
    number_of_sites = 10
    analyze_random_walk_on_circle(number_of_sites)