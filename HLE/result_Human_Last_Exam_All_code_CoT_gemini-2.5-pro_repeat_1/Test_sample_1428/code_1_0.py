import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes the 1D random walk on a circle, computes the transition matrix,
    verifies its eigenvalues/eigenvectors, and finds the rate of relaxation.
    """
    # Set the number of sites on the circle. N must be > 2 for a non-trivial result.
    N = 5

    # --- Step 1: Explain the Model and One-Step Transformation ---
    print("1. One-Dimensional Random Walk on a Circle (N sites)")
    print("------------------------------------------------------")
    print(f"We consider a symmetric random walk on a circle with N = {N} sites.")
    print("The one-step transformation of the probability distribution P_t is:")
    print("P_t+1(i) = 0.5 * P_t(i-1) + 0.5 * P_t(i+1)  (indices are modulo N)\n")

    # --- Step 2: Construct and Print the Transition Matrix A ---
    print(f"2. Transition Probability Matrix A (for N={N})")
    print("------------------------------------------------------")
    print("A_ij is the probability of moving from site j to site i in one step.")
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5
    print("The matrix A is:")
    print(A)
    print("\n")

    # --- Step 3: Verify the Eigenvectors and Eigenvalues ---
    print("3. Eigenvectors and Eigenvalues Verification")
    print("------------------------------------------------------")
    print("The eigenvectors v_n have components v_n[j] = exp(i * k_n * j),")
    print("with k_n = 2*pi*n/N. The eigenvalues are lambda_n = cos(k_n).\n")
    print(f"Let's verify this for n=1 and N={N}:")

    n = 1
    k_n = 2 * np.pi * n / N
    lambda_n_analytical = np.cos(k_n)
    v_n = np.exp(1j * k_n * np.arange(N))

    # The equation to verify is A * v_n = lambda_n * v_n
    print(f"For n = {n}:")
    print(f"k_{n} = (2 * pi * {n}) / {N} = {k_n:.5f}")
    print(f"Predicted eigenvalue: lambda_{n} = cos({k_n:.5f}) = {lambda_n_analytical:.5f}\n")

    # Calculate Left-Hand Side (LHS) and Right-Hand Side (RHS) of the equation
    lhs = A @ v_n
    rhs = lambda_n_analytical * v_n

    print("We check the equation: A * v_n = lambda_n * v_n")
    print("LHS (A * v_n) result:")
    print(np.round(lhs, 5))
    print("\nRHS (lambda_n * v_n) result:")
    print(np.round(rhs, 5))

    # Check if they are numerically close
    if np.allclose(lhs, rhs):
        print("\nVerification successful: The vector is an eigenvector with the predicted eigenvalue.\n")
    else:
        print("\nVerification failed.\n")

    # --- Step 4: Find the Rate of Relaxation ---
    print("4. Rate of Relaxation")
    print("------------------------------------------------------")
    print("The rate of relaxation to the stationary distribution is determined by")
    print("the second-largest eigenvalue of A, as it corresponds to the slowest decaying mode.\n")

    # The largest eigenvalue is always lambda_0 = cos(0) = 1.
    # The second-largest is lambda_1.
    lambda_second_largest_analytical = np.cos(2 * np.pi / N)

    print("The analytical formula for the second-largest eigenvalue is lambda_1 = cos(2*pi/N).")
    print("The final equation for the rate is:")
    print(f"Rate = cos(2 * pi / N) = cos(2 * {np.pi:.5f} / {N}) = cos({2*np.pi/N:.5f}) = {lambda_second_largest_analytical:.5f}\n")

    # For confirmation, we compute eigenvalues numerically from matrix A
    eigenvalues_numerical = np.linalg.eigvals(A)
    # Sort real parts in descending order to find the second largest
    eigenvalues_sorted = np.sort(np.real(eigenvalues_numerical))[::-1]
    lambda_second_largest_numerical = eigenvalues_sorted[1]

    print("For comparison, the numerically computed second-largest eigenvalue is:")
    print(f"{lambda_second_largest_numerical:.5f}")

if __name__ == '__main__':
    solve_random_walk_on_circle()