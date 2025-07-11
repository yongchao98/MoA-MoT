import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes a 1D random walk on a circle with N sites.
    """
    # Set the number of sites on the circle. You can change this value.
    N = 10

    # --- 1. One-step transformation of the probability distribution ---
    print("Step 1: One-step transformation")
    print("Let P_i(t) be the probability of being at site i at time t.")
    print("From any site, the walker moves to one of its two neighbors with probability 0.5.")
    print("The one-step transformation is given by the equation:")
    print("P_i(t+1) = 0.5 * P_{i-1}(t) + 0.5 * P_{i+1}(t)")
    print("where indices are cyclic (e.g., for i=1, i-1 is N; for i=N, i+1 is 1).")
    print("In vector form, p(t+1) = A * p(t), where A is the transition matrix.")
    print("-" * 50)

    # --- 2. Construct the Transition Probability Matrix A ---
    print("Step 2: Transition Probability Matrix A")
    # Initialize an N x N zero matrix
    # The element A[i, j] represents the probability of transitioning FROM j TO i.
    # The equation P_i = 0.5*P_{i-1} + 0.5*P_{i+1} means that the probability at state i
    # is sourced from states i-1 and i+1.
    # So, the non-zero elements in row i are A[i, i-1] = 0.5 and A[i, i+1] = 0.5.
    A = np.zeros((N, N))
    for i in range(N):
        # Contribution from the site to the "left" (index i-1)
        A[i, (i - 1 + N) % N] = 0.5
        # Contribution from the site to the "right" (index i+1)
        A[i, (i + 1) % N] = 0.5

    print(f"For N={N}, the transition matrix A is:")
    print(A)
    print("-" * 50)

    # --- 3. Verify Eigenvectors and Eigenvalues ---
    print("Step 3: Eigenvectors and Eigenvalues Verification")
    print(f"The theory states that the eigenvectors v_n have components (v_n)_j = exp(i * k_n * j)")
    print(f"where k_n = 2*pi*n/N, for n = 0, ..., N-1 and j = 0, ..., N-1.")
    print("The corresponding eigenvalue is lambda_n = cos(k_n).")

    # Let's verify this numerically for n=1
    n_test = 1
    k_test = 2 * np.pi * n_test / N
    lambda_test = np.cos(k_test)
    v_test = np.array([np.exp(1j * k_test * j) for j in range(N)])

    # Check if A @ v_test is equal to lambda_test * v_test
    lhs = A @ v_test
    rhs = lambda_test * v_test
    is_eigenvector = np.allclose(lhs, rhs)

    print(f"\nVerifying for n={n_test}:")
    print(f"Theoretical eigenvalue lambda_{n_test} = cos(2*pi/{N}) = {lambda_test:.6f}")
    print(f"Is A*v_{n_test} == lambda_{n_test}*v_{n_test}? -> {is_eigenvector}")

    # Compare all numerical eigenvalues with the theoretical formula
    numeric_eigenvalues, _ = np.linalg.eig(A)
    theoretical_eigenvalues = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    
    # Sort for comparison
    numeric_eigenvalues = np.sort(np.real(numeric_eigenvalues))
    theoretical_eigenvalues = np.sort(theoretical_eigenvalues)

    print(f"\nComparing all numerical and theoretical eigenvalues:")
    print(f"Verification successful: {np.allclose(numeric_eigenvalues, theoretical_eigenvalues)}")
    print("-" * 50)

    # --- 4. Rate of Relaxation ---
    print("Step 4: Rate of Relaxation")
    print("The system relaxes to its stationary distribution (corresponding to eigenvalue 1).")
    print("The rate of this relaxation is determined by the second-largest eigenvalue,")
    print("as it represents the slowest-decaying non-stationary mode.")

    # The eigenvalues are cos(2*pi*n/N).
    # The largest is for n=0: cos(0) = 1.
    # The second largest is for n=1 and n=N-1: cos(2*pi/N).
    lambda_second_largest_val = np.cos(2 * np.pi * 1 / N)

    print(f"\nThe largest eigenvalue is lambda_0 = cos(0) = 1.0")
    print("The second-largest eigenvalue determines the rate of relaxation.")
    
    print("\nThe final equation for the second-largest eigenvalue is:")
    n_second_largest = 1
    final_value = np.cos(2 * np.pi * n_second_largest / N)
    print(f"lambda_{n_second_largest} = cos(2 * pi * {n_second_largest} / {N}) = {final_value:.6f}")


solve_random_walk_on_circle()

<<<cos(2*pi/N)>>>