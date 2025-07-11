import numpy as np

def analyze_random_walk_on_circle(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function performs the following steps:
    1.  Describes the one-step transformation of the probability distribution.
    2.  Computes and prints the transition probability matrix A.
    3.  Analytically derives and numerically verifies the eigenvalues and eigenvectors of A.
    4.  Calculates the second-largest eigenvalue, which determines the rate of relaxation.

    Args:
        N (int): The number of sites on the circle.
    """
    if N <= 2:
        print("N must be greater than 2 for a non-trivial walk.")
        return

    print(f"Analyzing a 1D symmetric random walk on a circle with N = {N} sites.")
    
    # 1. One-step Transformation
    print("\n--- 1. One-step Transformation of the Probability Distribution ---")
    print("Let P_t(j) be the probability of being at site j (j=1,...,N) at time t.")
    print("For a symmetric random walk, the walker moves to an adjacent site with probability 1/2.")
    print("The one-step transformation is given by:")
    print("  P_{t+1}(j) = (1/2) * P_t(j-1) + (1/2) * P_t(j+1)")
    print("where the site indices are taken modulo N (e.g., for j=1, j-1 is N; for j=N, j+1 is 1).")
    print("In matrix form, this equation can be written as P_{t+1} = A * P_t, where A is the transition matrix.")

    # 2. Transition Probability Matrix A
    print("\n--- 2. The Transition Probability Matrix A ---")
    # We use 0-based indexing for the matrix (sites 0 to N-1)
    # A_ij is the probability of transitioning FROM j TO i.
    A = np.zeros((N, N))
    for j in range(N):
        # From site j, the walker can go to (j-1)%N or (j+1)%N
        A[(j - 1 + N) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5
    
    print("The transition matrix A is an NxN matrix where A_ij = P(i|j):")
    # Pretty print for smaller N
    if N <= 10:
        print(A)
    else:
        print(f"(Matrix is {N}x{N}, too large to display fully)")
        print(A[:4, :4])
        print("...")

    # 3. Eigenvectors and Eigenvalues
    print("\n--- 3. Eigenvectors and Eigenvalues ---")
    print("It can be shown that the eigenvectors v_n (for n=0, 1,...,N-1) have components:")
    print("  (v_n)_j = exp(i * j * k_n), where k_n = 2 * pi * n / N")
    print("\nTo show this, we apply the matrix A to the vector v_n:")
    print("  (A * v_n)_j = sum_l(A_jl * (v_n)_l)")
    print("              = (1/2)*(v_n)_{j-1} + (1/2)*(v_n)_{j+1}")
    print("              = (1/2)*exp(i*(j-1)*k_n) + (1/2)*exp(i*(j+1)*k_n)")
    print("              = exp(i*j*k_n) * [ (1/2)*exp(-i*k_n) + (1/2)*exp(i*k_n) ]")
    print("Using Euler's formula, cos(x) = (e^(ix) + e^(-ix))/2, the term in brackets is cos(k_n).")
    print("  (A * v_n)_j = exp(i*j*k_n) * cos(k_n) = (v_n)_j * cos(2*pi*n/N)")
    print("\nThis confirms that v_n is an eigenvector with a corresponding eigenvalue lambda_n:")
    print("  lambda_n = cos(2 * pi * n / N)")

    # 4. Rate of Relaxation and the Second-Largest Eigenvalue
    print("\n--- 4. Rate of Relaxation ---")
    print("The convergence to the stationary distribution is determined by the eigenvalues.")
    print("The largest eigenvalue is lambda_0 = cos(0) = 1, which corresponds to the non-decaying stationary state.")
    print("The rate of relaxation is determined by the magnitude of the next largest eigenvalue(s).")
    print("These are lambda_1 and lambda_{N-1}, since cos(x) is largest for x closest to 0.")
    print("  lambda_1 = lambda_{N-1} = cos(2 * pi / N)")
    
    # Calculation
    n = 1
    second_largest_eigenvalue = np.cos(2 * np.pi * n / N)

    print("\nThe equation for the second-largest eigenvalue is:")
    # "output each number in the final equation"
    print(f"  lambda_1 = cos( 2 * {np.pi:.5f} / {N} )")
    print(f"  lambda_1 = cos( {2 * np.pi / N:.5f} )")
    print(f"  lambda_1 = {second_largest_eigenvalue:.5f}")
    
    # Numerically compute eigenvalues from A to verify
    numerical_eigenvalues, _ = np.linalg.eig(A)
    # Sort in descending real part order
    sorted_numerical_eigenvalues = np.sort(numerical_eigenvalues.real)[::-1]
    
    print("\nVerification using numerical computation:")
    print(f"  Numerical second-largest eigenvalue: {sorted_numerical_eigenvalues[1]:.5f}")

    print("\nThe rate of relaxation is determined by this eigenvalue.")
    print(f"<<<{second_largest_eigenvalue}>>>")

# Execute the analysis for a circle of size N=10
if __name__ == '__main__':
    analyze_random_walk_on_circle(N=10)
