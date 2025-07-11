import sys
import math

def solve():
    """
    Reads a matrix from stdin, finds the non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read the size of the matrix
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)

        # The center coordinate for 1-based indexing
        center_coord = (n + 1) // 2

        k_val = 0
        k_row = -1
        k_col = -1

        # Find the location of the non-zero digit 'k'
        for r in range(1, n + 1):
            line = sys.stdin.readline()
            if not line:
                continue
            row_values = list(map(int, line.split()))
            for c in range(1, n + 1):
                # lists are 0-indexed, so we use c-1
                val = row_values[c - 1]
                if val != 0:
                    k_val = val
                    k_row = r
                    k_col = c
        
        # Calculate the minimum number of moves (Manhattan distance)
        row_moves = abs(k_row - center_coord)
        col_moves = abs(k_col - center_coord)
        total_moves = row_moves + col_moves

        # The prompt asks to output the numbers in the "final equation".
        # Based on the specified output format "k r c z", we will print these four integers.
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential errors with input reading
        print(f"Error processing input: {e}", file=sys.stderr)

solve()