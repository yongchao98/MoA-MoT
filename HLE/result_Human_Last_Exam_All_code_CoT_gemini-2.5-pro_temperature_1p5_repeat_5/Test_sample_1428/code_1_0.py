import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes the random walk on a circle, builds its transition matrix,
    verifies its eigenvectors and eigenvalues, and finds the rate of relaxation.
    """
    # Step 1: Define the size of the circle (number of sites)
    N = 10
    print(f"Analyzing a 1D random walk on a circle with N = {N} sites.")
    print("-" * 50)

    # Step 2: Construct the transition matrix A
    # A[i, j] is the probability of moving from site j to site i.
    A = np.zeros((N, N))
    for j in range(N):
        # Probability 1/2 to move to the right neighbor (i = j+1)
        A[(j + 1) % N, j] = 0.5
        # Probability 1/2 to move to the left neighbor (i = j-1)
        A[(j - 1 + N) % N, j] = 0.5
    
    print("Transition Matrix A:")
    print(A)
    print("-" * 50)

    # Step 3: Verify the eigenvector and eigenvalue relationship
    # We will show that A * v_n = lambda_n * v_n for a specific n (e.g., n=1)
    n_verify = 1
    print(f"Verifying the eigenvector/eigenvalue for n = {n_verify}...")

    # The proposed eigenvector v_n has components v_n(j) = exp(i * j * k_n)
    j_indices = np.arange(N)
    k_n = 2 * np.pi * n_verify / N
    v_n = np.exp(1j * j_indices * k_n)

    # The corresponding eigenvalue is lambda_n = cos(k_n)
    lambda_n = np.cos(k_n)

    # Calculate LHS: A @ v_n
    lhs = A @ v_n
    # Calculate RHS: lambda_n * v_n
    rhs = lambda_n * v_n

    print(f"Vector A @ v_{n_verify}:")
    print(lhs)
    print(f"Vector lambda_{n_verify} * v_{n_verify}:")
    print(rhs)
    
    # Check if they are close (allowing for floating point inaccuracies)
    is_verified = np.allclose(lhs, rhs)
    print(f"\nVerification successful: {is_verified}")
    print("-" * 50)

    # Step 4: Find the rate of relaxation
    print("Finding the rate of relaxation determined by the second-largest eigenvalue.")

    # Calculate all theoretical eigenvalues: lambda_n = cos(2*pi*n/N)
    all_n = np.arange(N)
    all_lambdas = np.cos(2 * np.pi * all_n / N)
    
    # Find the unique eigenvalues and sort them in descending order
    unique_lambdas = np.unique(all_lambdas)
    sorted_lambdas = np.sort(unique_lambdas)[::-1] # Sort descending

    largest_lambda = sorted_lambdas[0]
    second_largest_lambda = sorted_lambdas[1]

    print(f"\nAll unique eigenvalues sorted: {[f'{val:.4f}' for val in sorted_lambdas]}")
    
    # Print the final equation and its numerical result
    print("\nThe largest eigenvalue corresponds to the stationary distribution.")
    print(f"Its equation is lambda_0 = cos(2 * pi * 0 / N) = cos(0) = {largest_lambda:.4f}")
    
    print("\nThe rate of relaxation is determined by the second-largest eigenvalue.")
    print(f"The equation for the second-largest eigenvalue is: lambda_1 = cos(2 * pi / N)")
    print(f"For N={N}, the value is: cos(2 * {np.pi:.5f} / {N}) = {second_largest_lambda:.5f}")


if __name__ == "__main__":
    solve_random_walk_on_circle()