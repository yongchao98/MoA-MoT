import numpy as np

def solve_random_walk_on_circle(N):
    """
    Analyzes a 1D random walk on a circle with N sites.

    This function:
    1. Describes the one-step probability transformation.
    2. Computes and prints the transition probability matrix A.
    3. Calculates the rate of relaxation from the second-largest eigenvalue and prints the equation.
    """
    if not isinstance(N, int) or N <= 1:
        print("Error: N must be an integer greater than 1.")
        return

    # 1. The one-step transformation of the probability distribution p is p(t+1) = A * p(t).
    #    Let's compute the matrix A.
    print(f"--- Analysis for a circle with N = {N} sites ---")
    
    # Initialize an NxN zero matrix
    # A_ij is the probability of moving to state i from state j.
    A = np.zeros((N, N))

    # Populate the matrix based on the random walk rules.
    # The walker at site j moves to (j-1)%N or (j+1)%N with probability 0.5 each.
    for j in range(N):
        A[(j - 1) % N, j] = 0.5
        A[(j + 1) % N, j] = 0.5

    print("\n1. The one-step transformation of the probability distribution is p(t+1) = A * p(t),")
    print("where the transition probability matrix A is:")
    print(A)

    # 2. The rate of relaxation is determined by the second-largest eigenvalue.
    #    From the analytical solution, we know the eigenvalues are lambda_n = cos(2*pi*n/N).
    #    The largest eigenvalue is lambda_0 = 1.
    #    The second-largest eigenvalue is lambda_1 = cos(2*pi/N).

    rate_of_relaxation = np.cos(2 * np.pi / N)

    print("\n2. The rate of relaxation is determined by the second-largest eigenvalue, which is lambda_1.")
    print("\n   The calculation for this rate is as follows:")
    
    # Printing each number in the final equation as requested.
    print(f"   Rate = cos( 2 * pi / N )")
    print(f"   Rate = cos( 2 * {np.pi:.5f} / {N} )")
    print(f"   Rate = cos( {2 * np.pi / N:.5f} )")
    print(f"   Rate = {rate_of_relaxation:.5f}")


# --- Execution ---
# Set the number of sites on the circle. You can change this value.
number_of_sites = 10
solve_random_walk_on_circle(number_of_sites)
