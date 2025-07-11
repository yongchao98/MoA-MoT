import numpy as np
from functools import reduce
from operator import xor

def calculate_f_for_example(n, m):
    """
    Calculates the number of P-positions for an n x m matrix by brute force
    and determines the value of f(n, m).
    """
    p_pos_count = 0
    total_matrices = 1 << (n * m)

    # Iterate through all possible n x m binary matrices
    for i in range(total_matrices):
        # Create the matrix from the integer i
        matrix = np.zeros((n, m), dtype=int)
        temp_i = i
        for r in range(n):
            for c in range(m):
                if (temp_i & 1) == 1:
                    matrix[r, c] = 1
                temp_i >>= 1
        
        # A position is a P-position if the nim-sums of row and column counts are both 0.
        row_sums = np.sum(matrix, axis=1)
        col_sums = np.sum(matrix, axis=0)

        # The reduce function with an empty sequence would raise an error.
        # Handle cases where n or m is 0, though we assume n, m >= 1.
        nim_sum_rows = reduce(xor, row_sums) if n > 0 else 0
        nim_sum_cols = reduce(xor, col_sums) if m > 0 else 0

        if nim_sum_rows == 0 and nim_sum_cols == 0:
            p_pos_count += 1
            
    # The first player has a winning position with > 50% probability
    # if the number of P-positions is less than half the total.
    threshold = total_matrices / 2.0
    
    print(f"For a {n}x{m} matrix:")
    print(f"Total possible matrices = 2^({n}*{m}) = {total_matrices}")
    
    # This is the "final equation" we are evaluating
    print(f"Number of P-positions (losing for 1st player) = {p_pos_count}")
    print(f"Winning threshold (Total/2) = {threshold}")
    
    has_winning_majority = p_pos_count < threshold
    
    print(f"Is {p_pos_count} < {threshold}? {has_winning_majority}")
    
    f_value = 1 if has_winning_majority else 0
    print(f"Therefore, f({n}, {m}) = {f_value}")

# Run the demonstration for n=2, m=3
calculate_f_for_example(2, 3)