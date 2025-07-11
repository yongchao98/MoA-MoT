import numpy as np

def solve_random_walk_on_circle(N):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    if not isinstance(N, int) or N <= 2:
        print("Error: N must be an integer greater than 2.")
        return

    print(f"--- 1D Random Walk on a Circle with N = {N} sites ---\n")

    # 1. One-step transformation of the probability distribution
    print("1. One-Step Transformation:")
    print("The probability distribution evolves according to P_{t+1} = A * P_t, where:")
    print("P_{t+1}(j) = 0.5 * P_t((j-1) mod N) + 0.5 * P_t((j+1) mod N)\n")

    # 2. Compute the transition probability matrix A
    # We use indices 0 to N-1 for convenience with modulo arithmetic
    A = np.zeros((N, N))
    for j in range(N):
        A[(j - 1) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5

    print("2. The Transition Probability Matrix A:")
    # Pretty print the matrix if it's not too large
    if N <= 10:
        print(A)
    else:
        print(f"Matrix A is a {N}x{N} circulant matrix (too large to display).")
    print("\n")

    # 3. Verify the eigenvectors and eigenvalues
    print("3. Eigenvector and Eigenvalue Verification:")
    print("We verify that vectors v_n with components (v_n)_j = exp(i*j*k_n)")
    print("are eigenvectors with eigenvalues lambda_n = cos(k_n), where k_n = 2*pi*n/N.\n")

    all_verified = True
    theoretical_eigenvalues = []
    for n in range(N):
        k_n = 2 * np.pi * n / N
        lambda_n = np.cos(k_n)
        theoretical_eigenvalues.append(lambda_n)
        
        # Construct the eigenvector v_n
        v_n = np.array([np.exp(1j * j * k_n) for j in range(N)])
        
        # Check if A * v_n = lambda_n * v_n
        lhs = A @ v_n
        rhs = lambda_n * v_n
        if not np.allclose(lhs, rhs):
            all_verified = False
            print(f"Verification FAILED for n = {n}")
            break

    if all_verified:
        print("Successfully verified that v_n are eigenvectors for all n = 0, ..., N-1.")
        print("The corresponding eigenvalues are lambda_n = cos(2*pi*n/N).\n")

    # 4. Find the rate of relaxation
    print("4. Rate of Relaxation:")
    print("The rate of relaxation is determined by the second-largest eigenvalue.")
    
    # The largest eigenvalue is lambda_0 = 1.
    # The second-largest is lambda_1 = cos(2*pi/N).
    second_largest_eigenvalue = np.cos(2 * np.pi / N)

    print("The final equation for the second-largest eigenvalue is:")
    print(f"lambda_1 = cos(2 * pi / {N})")
    
    # As requested, outputting each number in the final equation
    print("\nComponent values in the equation:")
    print(f"Final value of lambda_1: {second_largest_eigenvalue}")
    print(f"Value of the integer 2: 2")
    print(f"Value of pi (approx): {np.pi}")
    print(f"Value of N: {N}")
    
    return second_largest_eigenvalue


if __name__ == '__main__':
    # Set the number of sites on the circle. You can change this value.
    N = 10
    rate = solve_random_walk_on_circle(N)
    # The final answer format requires the numerical value.
    if rate is not None:
        print(f"\n<<<__{rate}__>>>")
