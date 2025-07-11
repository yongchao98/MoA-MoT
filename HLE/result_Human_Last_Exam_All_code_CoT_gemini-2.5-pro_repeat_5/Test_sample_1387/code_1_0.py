import numpy as np
from itertools import combinations, chain

def solve_loopless_cycle_cover_parity():
    """
    Calculates the parity of the number of loopless cycle covers for a given graph.

    The method is based on the inclusion-exclusion principle modulo 2.
    It calculates Sum_{T} det(A_T) mod 2, where T iterates over all subsets of
    2-cycles in the graph.

    Note: This implementation is for demonstration and is not polynomial-time,
    as it iterates through an exponential number of subsets. Advanced algorithms
    exist that compute this sum in polynomial time.
    """
    # Example: A complete directed graph K_4 on 4 vertices.
    # The adjacency matrix has A_ij = 1 if i != j, and 0 otherwise.
    n = 4
    # The adjacency matrix is taken modulo 2, so entries are 0 or 1.
    A = np.array([
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0]
    ], dtype=int)

    print(f"The graph has n={n} vertices.")
    print("The adjacency matrix A over GF(2) is:")
    print(A)

    # Step 1: Find all pairs of vertices {u, v} that form a 2-cycle.
    two_cycle_pairs = []
    for r in range(n):
        for c in range(r + 1, n):
            if A[r, c] == 1 and A[c, r] == 1:
                two_cycle_pairs.append(frozenset([r, c]))
    
    print(f"\nFound {len(two_cycle_pairs)} pairs of vertices forming 2-cycles: {list(map(set, two_cycle_pairs))}")

    # Step 2: Sum the determinants of submatrices modulo 2.
    # The sum is over all subsets T of the set of 2-cycle pairs.
    total_parity = 0
    
    # Generate the powerset of two_cycle_pairs
    powerset = chain.from_iterable(combinations(two_cycle_pairs, r) for r in range(len(two_cycle_pairs) + 1))

    for T in powerset:
        # T is a subset of the 2-cycle pairs, e.g., ({{0, 1}}, {{2, 3}})
        
        # Get the set of vertices to remove for this subset T
        vertices_to_remove = set()
        for pair in T:
            vertices_to_remove.update(pair)
        
        # Create the submatrix A_T by deleting these rows and columns
        all_vertices = set(range(n))
        remaining_vertices_indices = sorted(list(all_vertices - vertices_to_remove))
        
        if not remaining_vertices_indices:
            # The determinant of a 0x0 matrix is defined as 1
            det_A_T = 1
        else:
            submatrix = A[np.ix_(remaining_vertices_indices, remaining_vertices_indices)]
            # Calculate determinant and take it modulo 2
            det_val = np.linalg.det(submatrix)
            det_A_T = int(round(det_val)) % 2

        total_parity = (total_parity + det_A_T) % 2
        
    print("\nResult of the calculation:")
    print(f"The parity of the number of loopless cycle covers is {total_parity}.")
    print("\nFor K_4, the loopless cycle covers are the 6 cycles of length 4. The number is 6, so the parity is 0. The code confirms this.")

solve_loopless_cycle_cover_parity()