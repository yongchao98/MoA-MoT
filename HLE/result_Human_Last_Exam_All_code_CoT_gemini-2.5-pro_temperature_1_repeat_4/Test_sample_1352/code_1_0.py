import numpy as np

def calculate_expected_R():
    """
    Calculates the expected number of rounds R for the switches to return to their initial state.
    """

    # Define influence sets (using 0-based indexing for convenience)
    influence_sets = {
        0: {1, 3, 5, 6},
        1: {2, 4, 5, 7},
        2: {3, 5},
        3: {4},
        4: {5, 7},
        5: {6},
        6: {7},
        7: {}
    }

    # Create the influence matrix A
    # A[j, i] = 1 if person i+1 influences person j+1
    A = np.zeros((8, 8), dtype=int)
    for person_idx, influenced_set in influence_sets.items():
        for influenced_idx in influenced_set:
            A[influenced_idx, person_idx] = 1

    def rank_F2(M):
        """
        Calculates the rank of a matrix over the field F_2 using Gaussian elimination.
        """
        mat = M.copy().astype(int)
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:
                    i += 1
                if i < rows:
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]
                    for k in range(rows):
                        if k != pivot_row and mat[k, j] == 1:
                            mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # Calculate powers of A
    A2 = (A @ A) % 2
    A4 = (A2 @ A2) % 2
    A8 = (A4 @ A4) % 2  # This should be the zero matrix

    # Matrices for which we need the null space size. M_p = T^p - I = A^p.
    matrices = {1: A, 2: A2, 4: A4, 8: A8}
    
    # Calculate N(p) = number of states with order dividing p
    num_vectors = {}
    for p, M_p in matrices.items():
        rank = rank_F2(M_p)
        dim_ker = 8 - rank
        num_vectors[p] = 2**dim_ker

    # Calculate C(p) = number of states with exact order p
    C = {}
    C[1] = num_vectors[1]
    C[2] = num_vectors[2] - num_vectors[1]
    C[4] = num_vectors[4] - num_vectors[2]
    C[8] = num_vectors[8] - num_vectors[4]

    # Calculate the sum of R over all initial states
    total_R_sum = C[1] * 1 + C[2] * 2 + C[4] * 4 + C[8] * 8
    
    # Calculate the expected value E[R]
    expected_R = total_R_sum / 256.0

    # Print the equation as requested
    print(f"E[R] = ({C[1]} * 1 + {C[2]} * 2 + {C[4]} * 4 + {C[8]} * 8) / 256 = {total_R_sum} / 256 ≈ {expected_R:.2f}")

calculate_expected_R()
<<<7.68>>>