import numpy as np

def solve_random_walk():
    """
    Analyzes the 1D random walk on a circle, derives its eigenvalues,
    and finds the rate of relaxation.
    """

    print("### 1. Markov Chain Model ###\n")
    print("Consider a 1D random walk on a circle with N sites.")
    print("Let P_j(t) be the probability of being at site j at time t.")
    print("At each step, the walker moves to an adjacent site (j-1 or j+1) with equal probability 1/2.")
    print("\nThe one-step transformation of the probability distribution is:")
    print("P_j(t+1) = (1/2) * P_{j-1}(t) + (1/2) * P_{j+1}(t)\n")
    print("This can be written in matrix form as p(t+1) = A * p(t), where p(t) is the vector of probabilities.")
    print("The transition matrix A is an N x N circulant matrix where A[j, j-1] = 1/2 and A[j, j+1] = 1/2 (indices mod N), and all other elements are 0.\n")

    print("### 2. Eigenvalue Derivation ###\n")
    print("We test if vectors v_n with components v_n(j) = exp(l*j*k_n) are eigenvectors,")
    print("where l is the imaginary unit and k_n = 2*pi*n/N for n = 0, 1, ..., N-1.\n")
    print("Let's compute the j-th component of the product A * v_n:")
    print("(A * v_n)_j = (1/2) * v_n(j-1) + (1/2) * v_n(j+1)")
    print("           = (1/2) * exp(l*(j-1)*k_n) + (1/2) * exp(l*(j+1)*k_n)")
    print("\nFactoring out exp(l*j*k_n):")
    print("(A * v_n)_j = exp(l*j*k_n) * (1/2) * [exp(-l*k_n) + exp(l*k_n)]")
    print("\nUsing Euler's formula, exp(ix) + exp(-ix) = 2*cos(x), the term in brackets simplifies:")
    print("(1/2) * [2*cos(k_n)] = cos(k_n)")
    print("\nSo, (A * v_n)_j = cos(k_n) * v_n(j).")
    print("This confirms that v_n are the eigenvectors, and the corresponding eigenvalues are:")
    print("λ_n = cos(2*pi*n / N)\n")

    print("### 3. Rate of Relaxation ###\n")
    print("The eigenvalues λ_n determine how the system evolves.")
    print("The largest eigenvalue is for n=0: λ_0 = cos(0) = 1. This corresponds to the stationary (uniform) distribution.")
    print("\nThe rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    print("This corresponds to n=1 (and n=N-1):")
    print("λ_1 = cos(2*pi*1 / N) = cos(2*pi / N)\n")
    
    print("### 4. Numerical Verification ###\n")
    N_example = 5
    print(f"Let's verify this for a specific case, N = {N_example}.\n")

    # Construct the transition matrix A
    A = np.zeros((N_example, N_example))
    for i in range(N_example):
        A[i, (i - 1) % N_example] = 0.5
        A[i, (i + 1) % N_example] = 0.5
    
    print("The transition matrix A for N=5 is:")
    print(A)

    # Calculate eigenvalues numerically
    eigenvalues = np.linalg.eigvals(A)
    # Sort them by magnitude in descending order
    eigenvalues_sorted = sorted(np.abs(eigenvalues), reverse=True)
    
    print("\nNumerically computed eigenvalues (magnitudes, sorted):")
    print(eigenvalues_sorted)

    # Calculate the theoretical second-largest eigenvalue
    lambda_1_theoretical = np.cos(2 * np.pi / N_example)
    print(f"\nTheoretical second-largest eigenvalue λ_1 = cos(2*pi/{N_example}) = {lambda_1_theoretical:.6f}")
    
    # In numpy's output, there might be two eigenvalues with this magnitude
    numerical_lambda_1 = eigenvalues_sorted[1]
    print(f"Numerical second-largest eigenvalue = {numerical_lambda_1:.6f}")
    
    if np.isclose(lambda_1_theoretical, numerical_lambda_1):
        print("\nThe theoretical and numerical results match.\n")
    else:
        print("\nDiscrepancy found between theoretical and numerical results.\n")
        
    print("---------------------------------------------------------------")
    print("### Final Answer ###")
    print("The rate of relaxation is determined by the second-largest eigenvalue.")
    print("The final equation for this eigenvalue is:")
    # The prompt requests to output each number in the final equation.
    print("lambda_1 = cos( (2 * pi) / N )")
    print("---------------------------------------------------------------")


if __name__ == '__main__':
    solve_random_walk()