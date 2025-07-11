import numpy as np

def solve_random_walk_on_circle():
    """
    Analyzes a 1D random walk on a circle, calculates its transition matrix,
    verifies its eigenvalues, and finds the rate of relaxation.
    """
    # Get N from user
    try:
        N_str = input("Enter the number of sites (N) on the circle (an integer > 1): ")
        N = int(N_str)
        if N <= 1:
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter an integer greater than 1. Using N=10 as default.")
        N = 10

    # Step 1: One-step transformation (for explanation)
    print("\n--- 1. One-Step Transformation ---")
    print("The probability distribution pi_t over the sites transforms in one step as:")
    print(f"pi_{{t+1}}(j) = 0.5 * pi_t(j-1) + 0.5 * pi_t(j+1) (indices are modulo {N})")

    # Step 2: Transition Matrix A
    print(f"\n--- 2. Transition Matrix A for N={N} ---")
    # A[i, j] is the probability of transition from j to i.
    A = np.zeros((N, N))
    for j in range(N):
        A[(j - 1) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5
    print("The transition matrix A is:")
    print(A)

    # Step 3: Eigenvalues
    print("\n--- 3. Eigenvalues Verification ---")
    print("The eigenvalues are given by the formula: lambda_n = cos(2 * pi * n / N)")
    # We can verify this numerically. Let's find the eigenvalues of the matrix A we constructed.
    numeric_eigenvalues, _ = np.linalg.eig(A)
    # Sort in descending order for comparison
    numeric_eigenvalues = np.sort(np.real(numeric_eigenvalues))[::-1]

    # Calculate theoretical eigenvalues
    theoretical_eigenvalues = np.array([np.cos(2 * np.pi * n / N) for n in range(N)])
    # Sort in descending order for comparison
    theoretical_eigenvalues = np.sort(theoretical_eigenvalues)[::-1]

    print("\nTheoretical eigenvalues (sorted):")
    print(theoretical_eigenvalues)
    print("\nNumerically computed eigenvalues from matrix A (sorted):")
    print(numeric_eigenvalues)
    print(f"\nVerification: The numerically computed eigenvalues match the theoretical formula ({np.allclose(numeric_eigenvalues, theoretical_eigenvalues)}).")

    # Step 4: Rate of Relaxation
    print("\n--- 4. Rate of Relaxation ---")
    print("The rate of relaxation to the stationary distribution is determined by the second-largest eigenvalue.")
    print("The largest eigenvalue is always lambda_0 = cos(0) = 1.")
    print("The second-largest eigenvalue corresponds to n=1, which is lambda_1 = cos(2 * pi / N).")

    # Final calculation as requested
    n = 1
    pi_val = np.pi
    arg_val = 2 * pi_val * n / N
    result = np.cos(arg_val)

    print("\nFinal calculation for the second-largest eigenvalue:")
    print(f"Equation: lambda_1 = cos(2 * pi * n / N)")
    print(f"Substituting n = {n} and N = {N}:")
    print(f"lambda_1 = cos(2 * {pi_val:.6f} * {n} / {N})")
    print(f"lambda_1 = cos({arg_val:.6f})")
    print(f"lambda_1 = {result:.6f}")

    print(f"\nThis eigenvalue, lambda_1, determines the rate of relaxation.")
    print(f"\n<<< {result} >>>")

if __name__ == '__main__':
    solve_random_walk_on_circle()