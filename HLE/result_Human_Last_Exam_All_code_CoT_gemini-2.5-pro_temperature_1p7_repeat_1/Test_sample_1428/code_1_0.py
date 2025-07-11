import numpy as np

def analyze_random_walk_on_circle(N=10):
    """
    Analyzes a 1D random walk on a circle with N sites.
    It prints the derivation of the transition matrix, its eigenvalues,
    and the relaxation rate determined by the second-largest eigenvalue.
    """
    
    print("--- Step 1: System Definition and Transition Matrix ---")
    print(f"We consider a 1D random walk on a circle with N = {N} sites.")
    print("At each step, the walker moves to one of its two neighbors with equal probability 1/2.")
    print("Let P_j be the probability of being at site j (j=0, 1, ..., N-1).")
    print("The one-step transformation of the probability distribution is:")
    print("P_j(t+1) = (1/2) * P_{j-1}(t) + (1/2) * P_{j+1}(t)  (indices are modulo N)")
    print("\nThis can be written in matrix form as P(t+1) = A * P(t), where A is the transition matrix.")

    # Construct the transition matrix A
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i - 1) % N] = 0.5
        A[i, (i + 1) % N] = 0.5

    print(f"\nThe transition probability matrix A for N={N} is:")
    # Set print options for better matrix display
    if N <= 20:
        print(A)
    else:
        print(f"(Matrix is {N}x{N}, too large to display)")

    print("\n\n--- Step 2: Eigenvectors and Eigenvalues Derivation ---")
    print("We show that vectors v_n with components (v_n)_j = e^(i * k_n * j) are eigenvectors of A,")
    print(f"where j = 0, ..., {N-1} and k_n = 2 * pi * n / N for n = 0, ..., {N-1}.")
    
    print("\nThe j-th component of the vector (A * v_n) is:")
    print("(A * v_n)_j = Sum_{l=0 to N-1} A_jl * (v_n)_l = (1/2) * ((v_n)_{j-1} + (v_n)_{j+1})")
    
    print("\nSubstituting the definition (v_n)_j = e^(i * k_n * j):")
    print("(A * v_n)_j = (1/2) * [e^(i * k_n * (j-1)) + e^(i * k_n * (j+1))]")
    print("Factoring out e^(i * k_n * j):")
    print("(A * v_n)_j = e^(i * k_n * j) * (1/2) * [e^(-i * k_n) + e^(i * k_n)]")

    print("\nUsing Euler's formula, e^(ix) + e^(-ix) = 2*cos(x), we get:")
    print("(A * v_n)_j = e^(i * k_n * j) * cos(k_n) = (v_n)_j * cos(k_n)")

    print("\nThis proves that A * v_n = lambda_n * v_n, where the eigenvalues are:")
    print(f"lambda_n = cos(k_n) = cos(2 * pi * n / N), for n = 0, 1, ..., {N-1}")

    print("\n\n--- Step 3: Numerical Verification for n=1 ---")
    # We will verify for n=1, a non-trivial case.
    n = 1
    k_n = 2 * np.pi * n / N
    print(f"For n = {n}, k_n = 2 * pi * {n} / {N} = {k_n:.4f} radians.")

    # Proposed eigenvector v_1
    j_indices = np.arange(N)
    v_n = np.exp(1j * k_n * j_indices)

    # Theoretical eigenvalue lambda_1
    lambda_n = np.cos(k_n)
    print(f"The theoretical eigenvalue is lambda_{n} = cos({k_n:.4f}) = {lambda_n:.4f}")

    # Numerically compute A * v_n and lambda_n * v_n to verify
    LHS = A @ v_n
    RHS = lambda_n * v_n

    # Check if they are close enough to be considered equal
    is_eigenvector = np.allclose(LHS, RHS)
    print(f"\nVerification: Is A * v_n numerically equal to lambda_n * v_n? {is_eigenvector}")
    assert is_eigenvector, "Numerical verification failed!"
    print("The vector is confirmed to be an eigenvector with the calculated eigenvalue.")

    print("\n\n--- Step 4: Rate of Relaxation ---")
    print("The rate of relaxation to the stationary distribution (which corresponds to the largest eigenvalue lambda_0 = 1) is determined by the second-largest eigenvalue.")
    print("The eigenvalues are lambda_n = cos(2 * pi * n / N).")
    
    print("\nTo find the second-largest eigenvalue, we look for n in {1, ..., N-1} that maximizes |cos(2*pi*n/N)|.")
    print("The cosine function is largest for arguments closest to 0 (or multiples of 2*pi).")
    print(f"The arguments are 2*pi*n/N. For n=1, the argument is 2*pi/N. For n=N-1, it is 2*pi*(N-1)/N = 2*pi - 2*pi/N.")
    print("cos(2*pi/N) and cos(2*pi - 2*pi/N) are equal and represent the largest value after cos(0)=1.")

    sl_eigenvalue_val = np.cos(2 * np.pi / N)
    
    print("\nTherefore, the second-largest eigenvalue, which determines the relaxation rate, is:")
    print("lambda_SL = cos(2 * pi / N)")
    
    print("\nFor this case, the equation with its numbers is:")
    final_equation = f"lambda_SL = cos(2 * pi / {N})"
    print(f"{final_equation} = cos(2 * {np.pi:.5f} / {N}) = {sl_eigenvalue_val:.8f}")

if __name__ == '__main__':
    # You can change the number of sites N here
    number_of_sites = 10
    analyze_random_walk_on_circle(N=number_of_sites)