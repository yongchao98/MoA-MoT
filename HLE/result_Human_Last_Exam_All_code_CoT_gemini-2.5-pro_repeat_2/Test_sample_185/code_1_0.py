import sys

def solve():
    """
    This function reads a matrix description from standard input,
    finds the non-zero element, and calculates the minimum moves
    to move it to the center.
    """
    try:
        # Read the size of the matrix
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())

        # The center coordinate for 1-based indexing
        center = n // 2 + 1

        k_val = 0
        k_row = -1
        k_col = -1

        # Read the matrix row by row to find the non-zero digit
        for r in range(1, n + 1):
            line = sys.stdin.readline().strip().split()
            for c in range(1, n + 1):
                val = int(line[c - 1])
                if val != 0:
                    k_val = val
                    k_row = r
                    k_col = c
                    # Found the unique non-zero digit, no need to search further
                    break
            if k_val != 0:
                break
        
        # Calculate the Manhattan distance to the center
        row_moves = abs(k_row - center)
        col_moves = abs(k_col - center)
        total_moves = row_moves + col_moves

        # Print the final result in the format: k r c z
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential input errors
        print(f"An error occurred: {e}", file=sys.stderr)

solve()
