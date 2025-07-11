import numpy as np

def solve_random_walk_on_circle(N):
    """
    Analyzes the random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    print(f"--- Analysis for a Random Walk on a Circle with N={N} sites ---")

    # --- 1. The Transformation and Transition Matrix ---
    print("\n[Part 1: The Transition Matrix A]")
    print("The probability of being at site j at time t+1 is given by:")
    print("P_{t+1}(j) = 0.5 * P_t(j-1) + 0.5 * P_t(j+1)")
    print("This corresponds to the transition matrix A where:")
    print("A[i, j] = 0.5 if i and j are neighbors, and 0 otherwise.\n")

    # Construct the transition matrix A
    A = np.zeros((N, N))
    for j in range(N):
        # Move right: j -> j+1
        A[(j + 1) % N, j] = 0.5
        # Move left: j -> j-1
        A[(j - 1) % N, j] = 0.5

    print(f"The constructed transition matrix A for N={N} is:")
    print(A)

    # --- 2. Eigenvectors and Eigenvalues Verification ---
    print("\n[Part 2: Eigenvector and Eigenvalue Verification]")
    print("The eigenvectors v_n have components (v_n)_j = exp(i * 2*pi*n*j / N)")
    print("The eigenvalues lambda_n are given by cos(2*pi*n / N)")
    
    # Let's verify for n=1
    n = 1
    print(f"\nVerifying for mode n={n}...")
    
    # Theoretical eigenvector v_1
    j_indices = np.arange(N)
    k_n = 2 * np.pi * n / N
    v_n = np.exp(1j * k_n * j_indices)

    # Theoretical eigenvalue lambda_1
    lambda_n_theory = np.cos(k_n)
    
    # Calculate A * v_n
    A_v_n = A @ v_n
    
    # Calculate lambda_n * v_n
    lambda_v_n = lambda_n_theory * v_n

    # Check if they are equal
    is_eigenvector = np.allclose(A_v_n, lambda_v_n)
    
    print(f"Theoretical eigenvalue lambda_{n}: {lambda_n_theory:.6f}")
    print(f"Check if A @ v_{n} equals lambda_{n} * v_{n}: {is_eigenvector}")
    if not is_eigenvector:
        print("Note: Discrepancy might be due to floating point inaccuracies.")

    # --- 3. Rate of Relaxation ---
    print("\n[Part 3: Rate of Relaxation]")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    
    # Numerically compute all eigenvalues from the matrix A
    eigenvalues_numeric = np.linalg.eigvals(A)
    # Sort them in descending order. They are real for a symmetric matrix.
    sorted_eigenvalues = np.sort(eigenvalues_numeric)[::-1]

    largest_eig_numeric = sorted_eigenvalues[0]
    second_largest_eig_numeric = sorted_eigenvalues[1] # There's a duplicate at index 2

    print(f"\nNumerically computed eigenvalues (sorted):")
    print(np.round(sorted_eigenvalues, 6))

    print(f"\nLargest eigenvalue (numerical): {largest_eig_numeric:.6f}")
    print("Theoretical largest eigenvalue (lambda_0 = cos(0)): 1.0")

    print(f"\nSecond-largest eigenvalue (numerical): {second_largest_eig_numeric:.6f}")
    
    # The final equation part as requested
    print("\nThe rate of relaxation is determined by the second-largest eigenvalue, whose equation is:")
    print("lambda_second_largest = cos( (2 * pi) / N )")
    
    # Compare numerical with theory
    second_largest_eig_theory = np.cos(2 * np.pi / N)
    print(f"\nValue from formula for N={N}: cos(2*pi/{N}) = {second_largest_eig_theory:.6f}")
    print(f"The numerical result matches the theoretical formula: {np.isclose(second_largest_eig_numeric, second_largest_eig_theory)}")
    print("\n" + "="*50)


# --- Main execution ---
# You can change N to see the results for different circle sizes.
N_sites = 5
solve_random_walk_on_circle(N_sites)