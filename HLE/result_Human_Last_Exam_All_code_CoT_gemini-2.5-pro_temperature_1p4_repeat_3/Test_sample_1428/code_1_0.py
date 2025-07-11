import numpy as np

def solve_random_walk_on_circle(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.

    Args:
        N (int): The number of sites on the circle.
    """
    if N <= 2:
        print("N must be greater than 2 for a non-trivial circle walk.")
        return

    # 1. One-step transformation description
    print(f"--- 1D Random Walk on a Circle with N={N} sites ---")
    print("\nThe probability P_i of being at site i evolves according to:")
    print("P_i(t+1) = 0.5 * P_{i-1}(t) + 0.5 * P_{i+1}(t)")
    print("(Indices are taken modulo N, using 0-based indexing from 0 to N-1)\n")

    # 2. Construct and print the transition matrix A
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5

    print("--- Transition Probability Matrix A ---")
    print(A)
    print("\n")

    # 3. Verify eigenvectors and eigenvalues
    print("--- Eigenvector and Eigenvalue Verification ---")
    print("We verify that v_n with components (v_n)_j = exp(i*j*k_n) are eigenvectors,")
    print("where k_n = 2*pi*n/N and the eigenvalue is lambda_n = cos(k_n).\n")

    all_verified = True
    for n in range(N):
        k_n = 2 * np.pi * n / N
        lambda_n_theory = np.cos(k_n)
        
        # Define the eigenvector v_n using 0-based indexing for sites j=0..N-1
        v_n = np.array([np.exp(1j * j * k_n) for j in range(N)])
        
        # Apply the transition matrix to the eigenvector
        Av = A @ v_n
        # Calculate the expected result
        lambda_v = lambda_n_theory * v_n
        
        # Check if A*v == lambda*v
        is_eigenvector = np.allclose(Av, lambda_v)
        if not is_eigenvector:
            all_verified = False
            
        print(f"For n = {n}:")
        print(f"  k_{n} = 2*pi*{n}/{N:.1f} = {k_n:.4f} rad")
        print(f"  Theoretical eigenvalue lambda_{n} = cos({k_n:.4f}) = {lambda_n_theory:.4f}")
        print(f"  Verification (A @ v_{n} == lambda_{n} * v_{n}): {is_eigenvector}")

    if all_verified:
        print("\nSuccessfully verified the analytical eigenvectors and eigenvalues for all n.\n")
    else:
        print("\nVerification failed for at least one n.\n")

    # 4. Rate of Relaxation
    print("--- Rate of Relaxation ---")
    # The largest eigenvalue corresponds to n=0
    lambda_max = np.cos(0)
    # The second-largest eigenvalue corresponds to n=1 (and n=N-1)
    lambda_second_largest = np.cos(2 * np.pi / N)

    print("The rate of relaxation to the stationary distribution is determined by the")
    print("second-largest eigenvalue of the transition matrix.")
    print("\nLargest eigenvalue (n=0):")
    print(f"lambda_0 = cos(2 * pi * 0 / {N}) = {lambda_max:.5f}")

    print("\nSecond-largest eigenvalue (n=1):")
    # The required output format: output each number in the final equation
    print(f"lambda_1 = cos(2 * pi / {N})")
    print(f"lambda_1 = {lambda_second_largest:.5f}")

if __name__ == '__main__':
    # You can change the number of sites N here
    number_of_sites = 10
    solve_random_walk_on_circle(number_of_sites)
    # The final answer for N=10 is cos(2*pi/10)
    final_answer = np.cos(2 * np.pi / number_of_sites)
    # print(f"\n<<< {final_answer:.5f} >>>") # Suppressing this for final output format