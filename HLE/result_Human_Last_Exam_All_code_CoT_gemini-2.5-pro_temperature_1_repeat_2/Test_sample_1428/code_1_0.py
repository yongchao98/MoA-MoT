import numpy as np

def analyze_random_walk_on_circle(N):
    """
    Analyzes the random walk on a circle with N sites.
    
    This function will:
    1. Construct the transition matrix A for the Markov chain.
    2. Define the eigenvector v_n and its corresponding eigenvalue lambda_n for n=1.
    3. Numerically verify the eigenvector equation A * v_n = lambda_n * v_n.
    4. Print the second-largest eigenvalue, which determines the relaxation rate.
    """
    if N < 3:
        print("N must be 3 or greater for a meaningful circle.")
        return None

    # Step 1: Construct the transition matrix A
    # A_ji is the probability of going from state i to state j.
    A = np.zeros((N, N))
    for i in range(N):
        # Move right: from i to (i+1)%N
        A[(i + 1) % N, i] = 0.5
        # Move left: from i to (i-1)%N
        A[(i - 1 + N) % N, i] = 0.5
    
    # Step 2: Define the eigenvector and eigenvalue for n=1
    # n=1 corresponds to the second-largest eigenvalue.
    n = 1
    # Wavenumber k_n = 2 * pi * n / N
    k_n = 2 * np.pi * n / N
    # Site indices j = [0, 1, ..., N-1]
    j = np.arange(N)
    # The eigenvector v_1 has components (v_1)_j = e^(i * k_n * j)
    v_n = np.exp(1j * k_n * j)
    
    # The corresponding eigenvalue lambda_1 = cos(k_n)
    lambda_n = np.cos(k_n)

    # Step 3: Verify the eigenvector equation A * v = lambda * v
    # This is the final equation we want to demonstrate.
    
    # Calculate the Left Hand Side (LHS) of the equation
    LHS = A @ v_n
    
    # Calculate the Right Hand Side (RHS) of the equation
    RHS = lambda_n * v_n
    
    # Step 4: Print all the components and the final result
    # Set print options for better readability
    np.set_printoptions(precision=4, suppress=True)
    
    print(f"Analysis for a random walk on a circle with N = {N} sites.")
    print("-" * 60)
    
    print("The one-step transition matrix A is:")
    print(A)
    print("-" * 60)

    print(f"We verify the eigenvector equation A * v = lambda * v for n={n}.")
    print(f"The eigenvector v_{n} (components e^(i*2*pi*{n}*j/{N})) is:")
    print(v_n)
    print("-" * 60)

    print(f"The corresponding eigenvalue lambda_{n} (cos(2*pi*{n}/{N})) is:")
    print(f"{lambda_n:.6f}")
    print("-" * 60)
    
    print(f"Let's compute each side of the equation A * v_{n} = lambda_{n} * v_{n}:")
    print(f"\nLHS (A * v_{n}):")
    print(LHS)
    
    print(f"\nRHS (lambda_{n} * v_{n}):")
    print(RHS)
    
    is_close = np.allclose(LHS, RHS)
    print(f"\nAre the LHS and RHS numerically equal? {is_close}")
    print("-" * 60)
    
    # Final Answer
    print("The rate of relaxation is determined by the second-largest eigenvalue.")
    print(f"The second-largest eigenvalue is lambda_1 = cos(2*pi/N).")
    
    return lambda_n

if __name__ == '__main__':
    # Set the number of sites on the circle
    N_sites = 10
    
    # Run the analysis
    second_largest_eigenvalue = analyze_random_walk_on_circle(N_sites)
    
    if second_largest_eigenvalue is not None:
        # Print the final numerical answer in the required format
        print(f"For N={N_sites}, the value is: {second_largest_eigenvalue:.6f}")
        print(f"<<<{second_largest_eigenvalue:.6f}>>>")
