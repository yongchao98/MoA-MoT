import sys
import math

def solve_matrix_problem():
    """
    This function reads a matrix from stdin, finds the non-zero element,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # 1. Read the matrix size N.
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)
        
        # Initialize variables to store the digit's value and 1-based coordinates.
        k_val = 0
        k_row = 0
        k_col = 0

        # 2. Read the matrix to find the non-zero digit and its coordinates.
        for i in range(n):
            # Read one row of the matrix as a list of strings.
            row_list = sys.stdin.readline().split()
            for j in range(n):
                if row_list[j] != '0':
                    k_val = int(row_list[j])
                    # Convert 0-based loop indices to 1-based coordinates.
                    k_row = i + 1
                    k_col = j + 1
                    
        # 3. Calculate the coordinate of the center.
        center = n // 2 + 1
        
        # 4. Calculate the minimum moves (Manhattan distance).
        row_moves = abs(k_row - center)
        col_moves = abs(k_col - center)
        total_moves = row_moves + col_moves
        
        # 5. Print the required integers: k, r, c, z.
        # The problem statement says "you still need to output each number in the final equation!".
        # The final equation is z = |r - center| + |c - center|.
        # The output below includes all the components: k, r, c and the final result z.
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError):
        # Gracefully handle empty input or invalid number format.
        pass

solve_matrix_problem()