import sys
import math

def solve_matrix_problem():
    """
    This function reads a matrix definition from standard input,
    finds the single non-zero digit, and calculates the minimum moves
    to move it to the center of the matrix.
    """
    try:
        # Read the size of the matrix, N.
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())

        k_value = 0
        start_row = -1
        start_col = -1

        # Read the N x N matrix row by row to find the non-zero digit.
        for i in range(n):
            line = sys.stdin.readline()
            if not line:
                continue
            
            # Split the row into numbers and find the non-zero one.
            row_values = list(map(int, line.strip().split()))
            for j in range(n):
                if row_values[j] != 0:
                    k_value = row_values[j]
                    # The problem uses 1-based indexing for rows and columns.
                    start_row = i + 1
                    start_col = j + 1
                    # Since there is only one non-zero digit, we could stop
                    # searching, but continuing to read all input is generally safer.

        # If no non-zero digit was found (which shouldn't happen based on the problem description)
        if start_row == -1:
            return

        # The center of an N x N matrix (with N being odd) is at (N/2 + 1, N/2 + 1).
        # We use integer division //.
        center_pos = n // 2 + 1

        # Calculate the minimum number of moves using the Manhattan distance.
        # z = |start_row - center_pos| + |start_col - center_pos|
        row_moves = abs(start_row - center_pos)
        col_moves = abs(start_col - center_pos)
        total_moves = row_moves + col_moves

        # Print the final result: k r c z
        # where r and c are the original coordinates and z is the total moves.
        print(f"{k_value} {start_row} {start_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential input errors gracefully.
        pass

solve_matrix_problem()
<<<5>>>