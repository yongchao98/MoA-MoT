import math
import random

def np_oracle_is_not_rigid(matrix, N, delta, r):
    """
    A dummy NP oracle.
    This function decides if a matrix M is NOT (delta, r)-rigid.
    This problem is in NP. A real implementation would use a SAT solver.
    The verifier would guess a low-rank matrix M' and check if the
    Hamming distance to M is small.
    
    For this demonstration, to ensure the main loop terminates,
    we'll simulate the fact that a random matrix is very likely to be rigid.
    We'll make this oracle almost always return False.
    """
    # A random matrix is rigid with very high probability.
    # We simulate this by having the oracle return True (not rigid) with a tiny probability.
    if random.random() < 0.001:
        # This matrix happens to be non-rigid
        return True
    else:
        # This matrix is rigid
        return False

def construct_rigid_matrix(N, delta):
    """
    Constructs a rigid matrix using a randomized algorithm with an NP oracle.
    """
    print(f"Searching for an {N}x{N} rigid matrix...")
    
    # Calculate the rank 'r' for which rigidity is expected.
    # The gap N-r is Omega(sqrt(N)). We use 2*sqrt(N) as a concrete example.
    if N < 4:
        # The formula for r might not make sense for very small N
        print("N is too small for a meaningful demonstration.")
        return None, 0
        
    r_float = N - 2 * math.sqrt(N)
    r = math.floor(r_float)
    
    # The number of entries that can be changed is at most k.
    k = math.floor(delta * N**2)

    count = 0
    while True:
        count += 1
        # 1. Generate a random N x N matrix over GF(2)
        m = [[random.randint(0, 1) for _ in range(N)] for _ in range(N)]
        
        # 2. Check if the matrix is rigid using the NP oracle
        # The oracle returns true if the matrix is NOT rigid.
        if not np_oracle_is_not_rigid(m, N, delta, r):
            print(f"Found a rigid matrix after {count} attempt(s).")
            # This matrix M is (delta, r)-rigid with high probability.
            return m, r
        else:
            print(f"Attempt {count}: Matrix was not rigid. Retrying...")

def main():
    """
    Main function to run the demonstration.
    """
    N = 100  # Size of the matrix
    delta = 0.1 # Fraction of entries that can be changed

    # The algorithm to find a rigid matrix
    matrix, r_val = construct_rigid_matrix(N, delta)

    if matrix:
        # This part of the code demonstrates the calculation of r
        # and prints the final equation as requested.
        print("\n" + "="*40)
        print("Rigidity Calculation:")
        print(f"For N = {N} and delta = {delta}, we can construct a (delta, r)-rigid matrix.")
        
        # Using the formula r = N - 2*sqrt(N)
        sqrt_n_val = math.sqrt(N)
        gap_val = 2 * sqrt_n_val
        final_r_val = N - gap_val
        
        print("The largest achievable r has the form: r = N - Theta(sqrt(N))")
        print("Using a concrete formula r = N - 2 * sqrt(N):")
        
        # Print the equation with each number.
        print("\nFinal Equation:")
        print(f"r = {N} - 2 * sqrt({N}) = {N} - 2 * {sqrt_n_val:.2f} = {N} - {gap_val:.2f} = {final_r_val:.2f}")
        print(f"So, we can construct a ({delta}, {math.floor(final_r_val)})-rigid matrix.")
        print("="*40)

if __name__ == "__main__":
    main()
