import numpy as np

def analyze_random_walk_on_circle(N=5):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function demonstrates the key properties of the random walk as a Markov chain:
    1. Defines the one-step transformation and the transition matrix A.
    2. Verifies the analytical eigenvectors and eigenvalues of A.
    3. Determines the rate of relaxation from the second-largest eigenvalue.

    Args:
        N (int): The number of sites on the circle.
    """
    if not isinstance(N, int) or N <= 2:
        print("Error: N must be an integer greater than 2.")
        return

    # Use 0-based indexing (0, 1, ..., N-1) for convenience in Python.
    # The physical results are independent of this choice.
    
    # --- 1. One-Step Transformation ---
    print(f"--- Analysis of 1D Random Walk on a Circle with N={N} sites ---")
    print("\nLet P_t(j) be the probability of being at site j at time t.")
    print(f"The sites are indexed j = 0, 1, ..., {N-1}.")
    print("For a symmetric random walk, the probability of moving to an adjacent site is 1/2.")
    
    print("\n[Step 1] One-Step Transformation of Probability Distribution:")
    print("The probability distribution evolves according to the equation:")
    print("P_{t+1}(j) = 0.5 * P_t((j-1) mod N) + 0.5 * P_t((j+1) mod N)")
    print("This can be written in matrix form as P_vector_{t+1} = A * P_vector_t.")

    # --- 2. Transition Matrix A ---
    print("\n[Step 2] Transition Probability Matrix A:")
    print("A[i, j] is the probability of transitioning from site j to site i.")
    # A[i, j] is non-zero only if j is a neighbor of i.
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5
    
    np.set_printoptions(precision=4, suppress=True)
    print(A)

    # --- 3. Eigenvectors and Eigenvalues Verification ---
    print("\n[Step 3] Eigenvectors and Eigenvalues Verification:")
    print("The analytical eigenvectors are v_n, with components v_n[j] = exp(i * j * 2*pi*n/N).")
    print("The corresponding analytical eigenvalues are lambda_n = cos(2*pi*n/N).")
    print("\nWe will now numerically verify the equation A * v_n = lambda_n * v_n for each n.")

    eigenvalues = []
    for n in range(N):
        print(f"\n{'='*15} Verifying for n = {n} {'='*15}")
        k_n = 2 * np.pi * n / N
        
        # Analytical eigenvalue
        lambda_n = np.cos(k_n)
        eigenvalues.append(lambda_n)
        
        # Proposed eigenvector
        v_n = np.array([np.exp(1j * j * k_n) for j in range(N)])
        
        # Calculate the left side of the equation: A * v_n
        Av_n = A @ v_n
        
        # Calculate the right side of the equation: lambda_n * v_n
        lambda_v_n = lambda_n * v_n
        
        print(f"Analytical Eigenvalue: lambda_{n} = cos(2*pi*{n}/{N}) = {lambda_n:.4f}")
        print("\n--- Equation Check: A * v_n = lambda_n * v_n ---")
        
        print("\nMatrix A:")
        print(A)
        print("\nEigenvector v_n:")
        print(v_n)
        print("\nResult of A * v_n (Left side):")
        print(Av_n)
        print(f"\nResult of lambda_{n} * v_n (Right side):")
        print(lambda_v_n)
        
        # Check if the two results are close enough
        if np.allclose(Av_n, lambda_v_n):
            print("\n[SUCCESS] The left side and right side are equal.")
        else:
            print("\n[FAILURE] The sides are not equal.")
        print(f"{'='*45}")

    # --- 4. Rate of Relaxation ---
    print("\n[Step 4] Rate of Relaxation:")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    
    # Sort eigenvalues in descending order
    sorted_eigenvalues = np.sort(np.array(eigenvalues))[::-1]
    
    print("\nCalculated eigenvalues in descending order:")
    print(sorted_eigenvalues)
    
    lambda_0 = sorted_eigenvalues[0]
    print(f"\nThe largest eigenvalue is lambda_0 = {lambda_0:.4f}, which is always 1.")
    print("This corresponds to the stationary (uniform) distribution.")
    
    # The second-largest eigenvalue corresponds to n=1 (and n=N-1)
    second_largest_eigenvalue = np.cos(2 * np.pi / N)
    
    print(f"\nThe second-largest eigenvalue is lambda_1 = cos(2*pi/N).")
    print("The deviation from the stationary state decays proportionally to (lambda_1)^t.")
    print("Thus, the rate of relaxation is determined by this value.")
    
    print("\n--- Final Answer ---")
    print("The rate of relaxation is determined by the second-largest eigenvalue, which is given by the formula:")
    print(f"lambda_1 = cos(2 * pi / N)")
    print(f"For N={N}, the value is: cos(2 * pi / {N}) = {second_largest_eigenvalue:.6f}")


if __name__ == '__main__':
    # You can change the number of sites (N) here
    number_of_sites = 5
    analyze_random_walk_on_circle(number_of_sites)