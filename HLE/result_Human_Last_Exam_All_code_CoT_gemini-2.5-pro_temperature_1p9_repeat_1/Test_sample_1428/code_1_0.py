import numpy as np

def solve_random_walk():
    """
    Analyzes the 1D random walk on a circle with N sites.
    
    This function:
    1. Sets the number of sites, N.
    2. Constructs the transition matrix A.
    3. Finds the second-largest eigenvalue, which determines the relaxation rate.
    4. Numerically verifies the eigenvalue equation for this mode.
    5. Prints the full analysis, including the final equation numbers.
    """
    N = 5
    print(f"--- Analysis for a 1D Random Walk on a Circle with N = {N} sites ---")

    # 1. Construct the transition matrix A
    A = np.zeros((N, N))
    for j in range(N):
        # Probability 1/2 to move to j+1 (or 0 if j=N-1)
        A[(j + 1) % N, j] = 0.5
        # Probability 1/2 to move to j-1 (or N-1 if j=0)
        A[(j - 1) % N, j] = 0.5

    print("\n1. The transition probability matrix A is:")
    print(A)

    # 2. Identify the eigenvalue (lambda_1) and eigenvector (v_1) for relaxation
    # The relaxation rate is determined by the second-largest eigenvalue, n=1.
    n = 1
    k_n = 2 * np.pi * n / N
    lambda_n = np.cos(k_n)
    
    # Construct the eigenvector v_n
    j_indices = np.arange(N)
    v_n = np.exp(1j * j_indices * k_n)

    print(f"\n2. The rate of relaxation is determined by the second-largest eigenvalue.")
    print(f"This corresponds to n={n}, giving the eigenvalue lambda_{n} = cos(2*pi*{n}/{N}).")
    print(f"   lambda_{n} = {lambda_n:.6f}")
    
    # 3. Verify the eigenvalue equation: A * v_n = lambda_n * v_n
    # We will print the numbers for this final equation.
    
    lhs = A @ v_n
    rhs = lambda_n * v_n
    
    print("\n3. Verifying the eigenvalue equation A * v_n = lambda_n * v_n:")
    
    print("\nVector A * v_n (Left-Hand Side):")
    for val in lhs:
        print(f"  {val.real: .6f} + {val.imag: .6f}j")

    print("\nVector lambda_n * v_n (Right-Hand Side):")
    for val in rhs:
        print(f"  {val.real: .6f} + {val.imag: .6f}j")

    # Check if they are numerically close
    are_close = np.allclose(lhs, rhs)
    print(f"\nAre the Left-Hand Side and Right-Hand Side vectors numerically equal? {are_close}")

    # 4. Final Conclusion
    print("\n--- Conclusion ---")
    print("The analysis is verified. The rate of relaxation towards the uniform distribution")
    print(f"is governed by the second-largest eigenvalue of the transition matrix, which is:")
    print(f"lambda_1 = cos(2*pi/N) = {lambda_n}")
    
    # The final answer requested by the user prompt
    return lambda_n

# Run the solver and get the final value
final_answer = solve_random_walk()
# The required output format is <<<value>>>
# print(f"\n<<<{final_answer}>>>")

solve_random_walk()