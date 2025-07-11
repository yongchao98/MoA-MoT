import math

def NP_ORACLE(partial_matrix, N, delta, r):
    """
    This is a placeholder for a call to an NP oracle (e.g., a SAT solver).
    It solves the following NP-complete problem:
    
    Problem: CAN_BECOME_NON_RIGID
    Input: A partially filled NxN matrix P.
    Question: Does there exist a completion of P, call it M, that is NOT
              (delta, r)-rigid?
              
              (i.e., does there exist a completion M, a matrix S with at most
              delta*N*N non-zero entries, and a matrix M_low of rank at most r
              such that M + S = M_low?)

    A real implementation would convert this to a SAT formula and call a solver.
    For this demonstration, we assume it exists and returns True or False.
    This function cannot be implemented efficiently in reality unless P=NP.
    """
    # In a real FNP algorithm, this function makes a call to an NP oracle.
    # Since we cannot implement an NP oracle, we'll simulate a possible scenario.
    # This simulation does not produce a truly rigid matrix.
    # The logic of the main algorithm, however, demonstrates the FNP construction.
    
    # A placeholder simulation:
    # A real oracle's result would depend on the complex properties of the partial_matrix.
    # We'll return False for a sparse matrix to guide the algorithm.
    num_filled = sum(1 for row in partial_matrix for entry in row if entry is not None)
    num_ones = sum(1 for row in partial_matrix for entry in row if entry == 1)
    if num_filled > 0 and (num_ones / num_filled) < 0.1:
        return False
    return True


def construct_rigid_matrix(N, delta, r):
    """
    Constructs a rigid matrix using a greedy search guided by an NP oracle.
    This is an FNP algorithm.
    """
    print(f"Attempting to construct a ({delta}, {r})-rigid {N}x{N} matrix.")
    
    # We need r = o(N) for this algorithm to be provably correct.
    # Let's check if the provided r satisfies a loose version of this.
    if r >= N / (math.log(N, 2) if N > 1 else 1):
        print("Warning: The rank r may be too large for this FNP construction to be guaranteed.")
        print("The theory requires r = o(N), e.g., r < N/log(N).\n")

    # Initialize a partial matrix with None values
    partial_matrix = [[None for _ in range(N)] for _ in range(N)]

    # Determine each entry of the matrix one by one
    for i in range(N):
        for j in range(N):
            # Try setting the (i, j) entry to 0
            partial_matrix[i][j] = 0
            
            # Ask the oracle if this partial matrix can be completed to a NON-rigid one.
            # The oracle takes a copy to not modify the current state.
            oracle_input = [row[:] for row in partial_matrix]
            can_become_non_rigid_with_0 = NP_ORACLE(oracle_input, N, delta, r)

            if not can_become_non_rigid_with_0:
                # If setting the bit to 0 guarantees that ANY completion is RIGID,
                # then we safely choose 0.
                # The entry partial_matrix[i][j] is already 0.
                pass
            else:
                # If setting to 0 MIGHT lead to a non-rigid matrix, our strategy
                # is to choose 1 instead. The theoretical guarantee ensures that
                # if a rigid matrix exists, this greedy choice preserves a path to one.
                partial_matrix[i][j] = 1
    
    return partial_matrix

if __name__ == '__main__':
    N = 8  # Size of the matrix
    delta = 0.1 # At least 10% of entries must be changed

    # For the largest r, we choose a function that is o(N).
    # For example, N / log(N). For N=8, log2(N)=3.
    # A very slowly growing function is log*(N) (iterated logarithm).
    # log*(8) = 3. So r could be floor(8/3) = 2.
    # Let's use r = floor(N / log2(N)) for this example.
    r = math.floor(N / (math.log(N, 2) if N > 1 else 1))

    # The largest r we can construct a matrix for is of the form o(N).
    # We propose r = floor(N / log*(N)) as the answer.
    # For N=8, log*(8) = 3, so r = floor(8/3) = 2.
    # For N=16, log*(16) = 4, so r = floor(16/4) = 4.
    # For N=65536, log*(65536) = 5, so r = floor(65536/5) = 13107.
    
    print("The largest rank 'r' for which an FNP algorithm can construct a")
    print(f"({delta}, r)-rigid matrix is of the form r = o(N). A concrete example for")
    print(f"the largest such r is r = floor(N / log*(N)), where log* is the iterated logarithm.\n")

    final_matrix = construct_rigid_matrix(N, delta, r)

    print("\nConstructed Matrix (via simulated FNP algorithm):")
    for row in final_matrix:
        # Printing each number in the final result
        print(" ".join(map(str, row)))
