import numpy as np

def solve_random_walk_on_circle(N):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function:
    1. Defines the one-step probability transformation.
    2. Constructs the transition matrix A.
    3. Shows that e^(i*k_n*j) are eigenvectors and finds the eigenvalues.
    4. Calculates the rate of relaxation from the second-largest eigenvalue.
    """
    if N < 3:
        print("N must be 3 or greater for a non-trivial circle.")
        return

    print(f"--- Analysis of a 1D Random Walk on a Circle with N = {N} sites ---")

    # --- 1. One-Step Transformation of Probability ---
    print("\n1. One-Step Transformation:")
    print("The probability distribution P(t) over the sites evolves in time according to:")
    print("P_i(t+1) = 0.5 * P_{i-1}(t) + 0.5 * P_{i+1}(t)")
    print("This can be written in matrix form as P(t+1) = A @ P(t), where A is the transition matrix.")

    # --- 2. Transition Probability Matrix (A) ---
    print("\n2. Transition Probability Matrix (A):")
    # The element A[i, j] is the probability of transitioning to state i FROM state j.
    # To land on i, you must have been at a neighbor of i.
    # This means A[i, j] is 0.5 if j is a neighbor of i.
    # We use 0-based indexing for sites 0, 1, ..., N-1.
    A = np.zeros((N, N))
    for i in range(N):
        # A walker at site (i-1)%N can move to i.
        A[i, (i - 1 + N) % N] = 0.5
        # A walker at site (i+1)%N can move to i.
        A[i, (i + 1) % N] = 0.5
    print("The transition matrix A is:")
    print(A)

    # --- 3. Eigenvectors and Eigenvalues ---
    print("\n3. Eigenvectors and Eigenvalues:")
    print("The eigenvectors v^(n) have components v_j = exp(i * k_n * j), with k_n = 2*pi*n/N.")
    print("The corresponding eigenvalues are lambda_n = cos(k_n).")
    
    # Let's verify this for n=1
    n_test = 1
    k_n = 2 * np.pi * n_test / N
    j_indices = np.arange(N) # Note: problem uses j=1..N, code uses j=0..N-1. The result is identical.
    v_n = np.exp(1j * k_n * j_indices)
    
    # Calculate A @ v_n
    A_v = A @ v_n
    
    # Calculate lambda_n * v_n
    lambda_n_calc = np.cos(k_n)
    lambda_v = lambda_n_calc * v_n

    print(f"\nVerification for n = {n_test}:")
    print(f"Is A @ v^(1) == lambda_1 * v^(1)? -> {np.allclose(A_v, lambda_v)}")
    
    print("\nEigenvalues are calculated as lambda_n = cos(2*pi*n/N):")
    eigenvalues = []
    for n in range(N):
        lam_n = np.cos(2 * np.pi * n / N)
        eigenvalues.append(lam_n)
        print(f"lambda_{n:<2} = cos(2*pi*{n}/{N}) = {lam_n:.5f}")

    # --- 4. Rate of Relaxation ---
    print("\n4. Rate of Relaxation:")
    # The largest eigenvalue corresponds to n=0.
    lambda_max = eigenvalues[0] # cos(0) = 1
    # The second-largest eigenvalue corresponds to n=1 (and n=N-1).
    lambda_second_max = eigenvalues[1] # cos(2*pi/N)

    # The rate of relaxation is determined by the spectral gap: 1 - lambda_second_max
    relaxation_rate = 1 - lambda_second_max
    
    print("The rate of convergence to the stationary distribution is determined by the second-largest eigenvalue.")
    print(f"Largest eigenvalue (lambda_0): {lambda_max}")
    print(f"Second-largest eigenvalue (lambda_1): {lambda_second_max:.5f}")
    
    print("\nThe rate of relaxation (spectral gap) is 1 - lambda_second_max.")
    print(f"Rate = 1 - cos(2*pi/{N})")
    print(f"     = 1 - {lambda_second_max:.5f}")
    print(f"     = {relaxation_rate:.5f}")

    # For large N, this rate is approximately 2*pi^2 / N^2
    if N > 10:
        approx_rate = 2 * (np.pi**2) / (N**2)
        print(f"\nFor large N, this is approximately 2*pi^2/N^2 = {approx_rate:.5f}")


if __name__ == '__main__':
    # You can change the number of sites N here
    number_of_sites = 10
    solve_random_walk_on_circle(number_of_sites)
