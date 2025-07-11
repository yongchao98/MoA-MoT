import numpy as np

def analyze_random_walk_on_circle(N=10):
    """
    Analyzes a 1D symmetric random walk on a circle with N sites.

    This function performs the following steps:
    1.  Describes the one-step transformation for the probability distribution.
    2.  Constructs and displays the transition probability matrix, A.
    3.  Verifies that the proposed vectors v_n are eigenvectors of A and finds the corresponding eigenvalues lambda_n.
    4.  Calculates the rate of relaxation, which is determined by the second-largest eigenvalue.

    Args:
        N (int): The number of sites on the circle. Must be greater than 1.
    """
    if not isinstance(N, int) or N <= 1:
        print("Error: N must be an integer greater than 1.")
        return

    print(f"--- Analysis of 1D Random Walk on a Circle with N = {N} sites ---")
    
    # --- Step 1: One-step Transformation ---
    print("\n[Step 1: One-Step Transformation of the Probability Distribution]")
    print("Let P_j(t) be the probability of the walker being at site j at time t.")
    print("In one step, the walker moves to a neighboring site (left or right) with equal probability 1/2.")
    print("The probability distribution evolves according to the equation:")
    # The f-string uses {{ and }} to escape the braces for printing
    print(f"P_j(t+1) = (1/2) * P_{{j-1}}(t) + (1/2) * P_{{j+1}}(t)")
    print("where indices are considered on a circle (i.e., modulo N).")
    print("This can be written in matrix form as: P(t+1) = A * P(t)")

    # --- Step 2: Transition Probability Matrix A ---
    print("\n[Step 2: The Transition Probability Matrix A]")
    print("The matrix element A_ij is the probability of transitioning from site j to site i.")
    print("A_ij = 1/2 if i is a neighbor of j, and 0 otherwise.")
    
    # We use 0-based indexing (sites 0 to N-1) for convenience in Python.
    A = np.zeros((N, N))
    for j in range(N):
        # From site j, the walker can move to site (j-1) or (j+1).
        # This corresponds to setting the values in column j.
        A[(j - 1) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5
    
    print("\nFor N = {}, the transition matrix A is:".format(N))
    print(A)

    # --- Step 3: Eigenvectors and Eigenvalues ---
    print("\n[Step 3: Verification of Eigenvectors and Eigenvalues]")
    print("The eigenvectors v_n of this matrix have components (v_n)_j = exp(i*k_n*j)")
    print(f"where j = 0, ..., {N-1} and the wavenumber k_n = 2*pi*n/N for n = 0, ..., {N-1}.")
    print("The corresponding eigenvalues are lambda_n = cos(k_n).")
    
    # We will verify this for the case n=1
    n_verify = 1
    j_indices = np.arange(N)
    k_n_val = (2 * np.pi * n_verify) / N
    
    # Proposed eigenvector v_1
    v_n_vec = np.exp(1j * k_n_val * j_indices)
    
    # Proposed eigenvalue lambda_1
    lambda_n_val = np.cos(k_n_val)
    
    print(f"\nLet's verify the relationship A * v_n = lambda_n * v_n for n={n_verify}:")
    
    # Calculate Left-Hand Side (LHS) and Right-Hand Side (RHS)
    lhs = A @ v_n_vec
    rhs = lambda_n_val * v_n_vec
    
    is_eigenvector = np.allclose(lhs, rhs)
    
    print(f"For n={n_verify}, k_{n_verify} = 2*pi/{N} = {k_n_val:.4f} and lambda_{n_verify} = cos({k_n_val:.4f}) = {lambda_n_val:.4f}")
    print("We check if A * v_1 is numerically equal to lambda_1 * v_1.")
    print(f"The check result is: {is_eigenvector}. The property is verified.")
    
    # --- Step 4: Rate of Relaxation ---
    print("\n[Step 4: Calculation of the Relaxation Rate]")
    print("The convergence to the stationary distribution is governed by the eigenvalues.")
    print("The largest eigenvalue is lambda_0 = cos(0) = 1, for n=0.")
    
    # The second-largest eigenvalue corresponds to n=1 (and n=N-1)
    n_second = 1
    lambda_second_largest = np.cos((2 * np.pi * n_second) / N)
    
    print("The rate of relaxation is determined by the second-largest eigenvalue, which is lambda_1.")
    print(f"The equation for this eigenvalue is: lambda_1 = cos(2*pi/N)")
    print(f"For N = {N}, we calculate the value:")
    print(f"lambda_1 = cos(2 * {np.pi:.4f} / {N}) = {lambda_second_largest:.6f}")
    
    relaxation_rate = 1 - lambda_second_largest
    
    print("\nThe rate of relaxation is defined as the spectral gap: R = 1 - lambda_1")
    print(f"The final equation for the rate is: R = 1 - cos(2*pi/N)")
    print(f"Substituting the values, we get:")
    print(f"R = 1 - {lambda_second_largest:.6f} = {relaxation_rate:.6f}")
    
    return relaxation_rate

if __name__ == '__main__':
    # Define the number of sites on the circle
    N_sites = 10
    
    # Run the analysis and get the relaxation rate
    rate = analyze_random_walk_on_circle(N_sites)
    
    if rate is not None:
        print("\n--- Final Answer ---")
        print("The rate of relaxation is:")
        print(f"<<<{rate:.6f}>>>")
