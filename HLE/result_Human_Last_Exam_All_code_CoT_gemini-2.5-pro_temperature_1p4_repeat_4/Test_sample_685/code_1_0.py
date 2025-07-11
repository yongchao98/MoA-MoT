import math
import functools
import operator

def solve_nim_probability(n, m):
    """
    Calculates the exact probability that a random n x m binary matrix
    represents a P-position (losing position) in 2D-Generalized NIM.

    This is done by iterating through all possible 2^(n*m) matrices.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """
    if n * m > 20:
        print("The specified dimensions (n={}, m={}) are too large for an exhaustive search.".format(n, m))
        print("The complexity is O(n*m * 2^(n*m)), which is computationally infeasible.")
        return

    total_matrices = 1 << (n * m)
    p_position_count = 0

    # Iterate through all possible matrices represented by an integer from 0 to 2^(n*m)-1
    for i in range(total_matrices):
        matrix = []
        temp_i = i
        # Construct the n x m matrix from the integer i
        for r in range(n):
            row = []
            for c in range(m):
                row.append(temp_i & 1)
                temp_i >>= 1
            matrix.append(row)

        # Calculate row sums
        row_sums = [sum(row) for row in matrix]

        # Calculate column sums
        col_sums = [0] * m
        for r in range(n):
            for c in range(m):
                col_sums[c] += matrix[r][c]

        # Calculate nim-sums
        row_nim_sum = functools.reduce(operator.xor, row_sums)
        col_nim_sum = functools.reduce(operator.xor, col_sums)

        # Check for P-position condition
        if row_nim_sum == 0 and col_nim_sum == 0:
            p_position_count += 1
            
    probability_p_pos = p_position_count / total_matrices
    is_first_player_win_probable = probability_p_pos < 0.5
    f_nm = 1 if is_first_player_win_probable else 0

    print(f"For a {n}x{m} matrix:")
    print(f"Total number of possible matrices: {total_matrices}")
    print(f"Number of matrices that are P-positions: {p_position_count}")
    print(f"Probability of a P-position: {p_position_count} / {total_matrices} = {probability_p_pos:.4f}")
    print(f"Probability of an N-position (first player winning): {1 - probability_p_pos:.4f}")
    print(f"The value of f({n}, {m}) is {f_nm} because the P(N-position) > 0.5 is {is_first_player_win_probable}.")

# Example for a 2x2 matrix
solve_nim_probability(2, 2)
print("-" * 20)
# Example for a 2x3 matrix
solve_nim_probability(2, 3)
