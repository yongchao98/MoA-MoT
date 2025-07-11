import numpy as np

def analyze_random_walk_on_circle(N):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    if N < 3:
        print("N must be 3 or greater for a non-trivial circle.")
        return

    print(f"--- Analysis for a random walk on a circle with N = {N} sites ---\n")

    # Step 1: Define the one-step transformation
    print("1. One-step transformation of the probability distribution p(t):")
    print("p_j(t+1) = (1/2) * p_{(j-1) mod N}(t) + (1/2) * p_{(j+1) mod N}(t)")
    print("This can be written as p(t+1) = A * p(t), where A is the transition matrix.\n")

    # Step 2: Construct the transition matrix A
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5
    print("2. The transition probability matrix A is:")
    print(A)
    print("\n")

    # Step 3: Verify the eigenvectors and eigenvalues
    print("3. Verification of Eigenvectors and Eigenvalues:")
    print("The theoretical eigenvectors are v_n with components (v_n)_j = exp(i * k_n * j)")
    print("where i is the imaginary unit and k_n = 2*pi*n / N for n = 0, 1, ..., N-1.")
    print("The corresponding theoretical eigenvalues are lambda_n = cos(k_n) = cos(2*pi*n / N).\n")
    
    # Let's verify for n=1
    n_verify = 1
    k_n = 2 * np.pi * n_verify / N
    lambda_n_theory = np.cos(k_n)
    v_n = np.exp(1j * k_n * np.arange(N))
    
    # Calculate A * v_n
    A_v_n = A @ v_n
    # Calculate lambda_n * v_n
    lambda_v_n = lambda_n_theory * v_n

    print(f"Verifying for n = {n_verify}:")
    print(f"Theoretical eigenvalue lambda_{n_verify} = cos(2*pi/{N}) = {lambda_n_theory:.4f}")
    # print(f"Eigenvector v_{n_verify} =\n{np.round(v_n, 4)}\n")
    # print(f"A @ v_{n_verify} =\n{np.round(A_v_n, 4)}\n")
    # print(f"lambda_{n_verify} * v_{n_verify} =\n{np.round(lambda_v_n, 4)}\n")
    # Check if they are numerically close
    is_eigenvector = np.allclose(A_v_n, lambda_v_n)
    print(f"Is A*v_{n_verify} == lambda_{n_verify}*v_{n_verify}? {is_eigenvector}\n")

    # Step 4: Find the rate of relaxation
    print("4. Rate of Relaxation:")
    # Calculate all theoretical eigenvalues
    eigenvalues_theory = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    eigenvalues_theory_sorted = np.sort(eigenvalues_theory)[::-1]
    
    print(f"The theoretical eigenvalues lambda_n for n=0,...,{N-1} are:")
    print(np.round(eigenvalues_theory, 4))
    
    largest_eigenvalue = eigenvalues_theory_sorted[0]
    second_largest_eigenvalue = eigenvalues_theory_sorted[1]
    
    print(f"\nThe largest eigenvalue is lambda_0 = {largest_eigenvalue:.4f}, which corresponds to the stationary distribution.")
    
    # Note on periodicity for even N
    if N % 2 == 0:
        min_eigenvalue = eigenvalues_theory_sorted[-1]
        print(f"For even N, the eigenvalue lambda_{N//2} = {min_eigenvalue:.4f}. Since its magnitude is 1,")
        print("the chain is periodic and does not converge to a single stationary distribution.")
        print("However, the rate of mixing or diffusion is still governed by the second-largest eigenvalue.")

    print(f"\nThe second-largest eigenvalue is lambda_1 = cos(2*pi/N) = {second_largest_eigenvalue:.4f}.")
    print("The rate of relaxation to the stationary distribution is determined by the spectral gap,")
    print("which measures how quickly non-stationary components of the distribution decay.")
    print("Spectral Gap = 1 - (second-largest eigenvalue value)")
    print(f"Rate of relaxation = 1 - cos(2*pi/N)")

# Example usage with N=6
if __name__ == '__main__':
    N_sites = 6
    analyze_random_walk_on_circle(N_sites)