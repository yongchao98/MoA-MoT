import sys
import math

def solve():
    """
    Reads a matrix, finds a non-zero digit, and calculates the moves
    to bring it to the center.
    """
    try:
        # Read the size of the matrix
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())
    except (ValueError, IndexError):
        # Handle invalid input for N
        return

    k_val = 0
    k_row = -1
    k_col = -1

    # Read the matrix to find the non-zero digit and its 1-based position
    for i in range(1, n + 1):
        try:
            line = sys.stdin.readline()
            if not line:
                break
            row_vals = list(map(int, line.split()))
            
            # Only search if we haven't found the digit yet
            if k_val == 0:
                for j in range(n):
                    if row_vals[j] != 0:
                        k_val = row_vals[j]
                        k_row = i
                        k_col = j + 1
                        # Once found, we can stop searching the inner loop
                        break
        except (ValueError, IndexError):
            # Skip malformed rows
            continue

    # The target row and column (1-based)
    target_pos = n // 2 + 1

    # The final equation for moves is the Manhattan distance:
    # moves = abs(k_row - target_pos) + abs(k_col - target_pos)
    row_moves = abs(k_row - target_pos)
    col_moves = abs(k_col - target_pos)
    total_moves = row_moves + col_moves

    # The numbers in the final equation are k_val, k_row, k_col, and total_moves.
    # We print them as requested by the output format.
    print(f"{k_val} {k_row} {k_col} {total_moves}")

solve()