import sys

def solve():
    """
    Reads a matrix, finds the non-zero digit, and calculates the minimum moves
    to move it to the center.
    """
    try:
        # Read the size of the matrix
        n = int(sys.stdin.readline())

        k_val = -1
        k_row = -1
        k_col = -1

        # Read the matrix row by row to find the non-zero digit
        for i in range(n):
            row_str = sys.stdin.readline().split()
            for j in range(n):
                val = int(row_str[j])
                if val != 0:
                    k_val = val
                    # Using 1-based indexing for rows and columns as per the problem
                    k_row = i + 1
                    k_col = j + 1
                    break
            if k_val != -1:
                # To be more efficient, we could stop reading the rest of the matrix,
                # but for N < 17, reading the whole input is fine.
                pass
        
        # The center of the matrix (using 1-based indexing)
        center_coord = (n // 2) + 1

        # Calculate the Manhattan distance
        row_moves = abs(k_row - center_coord)
        col_moves = abs(k_col - center_coord)
        total_moves = row_moves + col_moves
        
        # Print the final result: k, original_row, original_col, total_moves
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()