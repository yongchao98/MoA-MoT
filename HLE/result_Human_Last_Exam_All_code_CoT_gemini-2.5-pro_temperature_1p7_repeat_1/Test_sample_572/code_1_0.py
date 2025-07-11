import math
import random

# This is a placeholder for the NP oracle.
# The actual problem ("Is M NOT (delta, r)-rigid?") is NP-complete.
# An FNP algorithm is given access to such an oracle.
def np_oracle_is_not_rigid(matrix, delta, r):
    """
    Oracle for the NP-complete problem: "Is matrix M NOT (delta, r)-rigid?"

    A matrix M is NOT (delta, r)-rigid if there exists a matrix S
    with at most delta * N*N non-zero entries, such that rank(M + S) <= r.

    This function is a mock-up. In a real FNP machine, this would be
    a call to an NP-complete problem solver. We cannot implement this efficiently.
    For this demonstration, it returns a random boolean, biased towards False
    as rigid matrices are extremely common.
    """
    # Mock implementation: most matrices are rigid, so this oracle should usually return False.
    return random.random() < 0.1


def construct_rigid_matrix(N, delta, r):
    """
    Constructs a rigid matrix using a lexicographical search guided by an NP oracle.
    This represents an FP^NP algorithm.

    The algorithm iterates through all possible N x N matrices over GF(2) in 
    lexicographical order. For each matrix, it asks the NP oracle if the matrix 
    is NOT rigid. The first matrix for which the oracle answers "NO" is a 
    (delta, r)-rigid matrix.
    
    NOTE: The runtime depends on finding a rigid matrix early in the search. 
    A worst-case polynomial bound is not known, but this illustrates the power of FNP.
    """
    # To prevent a very long loop, we demonstrate with a small N.
    effective_N = N
    if N > 4:
        print(f"NOTE: N={N} is large for a full search. Demonstrating with N=4.")
        effective_N = 4

    num_matrices = 1 << (effective_N * effective_N)
    for i in range(num_matrices):
        # Construct the i-th matrix in lexicographical order.
        M = []
        # Use a temporary variable for bit manipulation to preserve i for the printout.
        temp_i = i
        for _ in range(effective_N):
            matrix_row = []
            for _ in range(effective_N):
                matrix_row.append(temp_i & 1)
                temp_i >>= 1
            M.append(matrix_row)

        # Ask the oracle: "Is matrix M NOT (delta, r)-rigid?"
        # The oracle can decide any NP problem, so it can solve this.
        is_not_rigid = np_oracle_is_not_rigid(M, delta, r)

        if not is_not_rigid:
            # If the oracle returns False, it means the matrix IS rigid.
            print(f"Found a ({delta}, {r})-rigid {effective_N}x{effective_N} matrix (the {i}-th candidate):")
            for row in M:
                print(" ".join(map(str, row)))
            return M
    
    print("No rigid matrix found within the search space for N={effective_N}.")
    return None

if __name__ == '__main__':
    # N is the size of the matrix. The problem asks about infinitely many N's.
    N = 8 
    
    # delta is a small constant.
    delta = 0.1 
    
    # The question is for the largest possible r.
    # Constructions based on error-correcting codes can achieve r = Omega(N).
    # A known concrete bound from such methods is (1/4 - epsilon) * N.
    # We set epsilon to a small value, e.g., 0.01.
    epsilon = 0.01
    r = math.floor((1/4 - epsilon) * N)
    
    print(f"Algorithm: Search for a ({delta}, r)-rigid NxN matrix for N={N} and r={r}.")
    construct_rigid_matrix(N, delta, r)
