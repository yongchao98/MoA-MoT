import numpy as np

def random_walk_analysis():
    """
    Analyzes a 1D random walk on a circle, calculates the transition matrix,
    its eigenvalues, and determines the rate of relaxation.
    """
    # You can change N to any integer value >= 3
    N = 10

    print("--- Analysis of a 1D Random Walk on a Circle ---")
    print(f"Number of sites N = {N}\n")

    # --- Step 1: Mathematical Formulation and Transition Matrix ---
    print("Step 1: The Transition Matrix A")
    print("The random walk is a Markov chain. The probability distribution P_t evolves as P_{t+1} = A * P_t.")
    print("A particle at site j moves to j-1 or j+1 with probability 1/2 each.")
    print("The transition matrix A (size NxN) has elements A_ij = P(next_state=i | current_state=j).")
    print("This results in a symmetric circulant matrix where non-zero elements are A[i, (i-1)%N] = 0.5 and A[i, (i+1)%N] = 0.5.")
    
    # Constructing the Transition Matrix A (using A_ij = P(i|j) convention)
    A = np.zeros((N, N))
    for j in range(N):
        # From state j, one can transition to state (j-1)%N or (j+1)%N
        A[(j - 1 + N) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5
    
    print("\nConstructed Transition Matrix A:")
    print(A)
    print("\n")

    # --- Step 2: Eigenvectors and Eigenvalues ---
    print("Step 2: Eigenvectors and Eigenvalues")
    print("The eigenvectors of this matrix have components v_n[j] = exp(i * k_n * j),")
    print(f"with wave number k_n = 2 * pi * n / N for n = 0, 1, ..., N-1.")
    print("The corresponding eigenvalues are lambda_n = cos(k_n) = cos(2 * pi * n / N).")
    
    # Verification with numerical calculation
    # For a symmetric matrix, eigenvalues are real.
    eigenvalues_numerical = np.linalg.eigvals(A)
    # Sort in descending order
    eigenvalues_numerical.sort()
    eigenvalues_numerical = eigenvalues_numerical[::-1]

    print("\nNumerically calculated eigenvalues (sorted):")
    print(np.round(eigenvalues_numerical, 6))

    theoretical_eigs = sorted([np.cos(2 * np.pi * n / N) for n in range(N)], reverse=True)
    print("\nTheoretical eigenvalues (sorted):")
    print(np.round(theoretical_eigs, 6))
    
    # Check if they match
    assert np.allclose(eigenvalues_numerical, theoretical_eigs), "Mismatch between numerical and theoretical eigenvalues!"
    print("\nNumerical and theoretical eigenvalues match.\n")

    # --- Step 3: Rate of Relaxation ---
    print("Step 3: Rate of Relaxation")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    print("The largest eigenvalue is lambda_0 = cos(0) = 1.")
    print("The second-largest eigenvalue value is lambda_1, which corresponds to n=1.")
    
    # Identify the second-largest eigenvalue
    second_largest_eig_val = eigenvalues_numerical[1]
    
    # The final equation as requested
    n_val = 1
    k_val = 2 * np.pi * n_val / N
    
    print("\n--- Final Equation for the Second-Largest Eigenvalue ---")
    print(f"The second-largest eigenvalue lambda_1 is calculated for n = {n_val}:")
    
    # Print each number in the equation
    print(f"lambda_1 = cos(2 * pi * n / N)")
    print(f"lambda_1 = cos(2 * {np.pi:.4f} * {n_val} / {N})")
    print(f"lambda_1 = cos({k_val:.4f})")
    print(f"lambda_1 = {np.cos(k_val):.8f}")
    
    print("\nThis value matches the numerically computed second-largest eigenvalue:")
    print(f"Numerical value = {second_largest_eig_val:.8f}")

    return np.cos(k_val)

# Run the analysis
final_answer = random_walk_analysis()
print(f"\n<<<The rate of relaxation is determined by the second-largest eigenvalue, which is {final_answer:.8f}>>>")