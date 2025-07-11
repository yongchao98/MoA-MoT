import sys

def solve():
    """
    Reads an N x N matrix, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    # Read all input lines from stdin.
    # The first line is N, the rest is the matrix.
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
        n = int(lines[0])
        matrix_lines = lines[1:]
    except (IOError, ValueError, IndexError):
        # Handle potential empty input or malformed N
        return

    k_val = 0
    k_row = -1
    k_col = -1

    # Find the non-zero digit 'k' and its 1-based coordinates 'r' and 'c'.
    # We can break early once the number is found.
    found = False
    for i in range(n):
        try:
            row_vals = list(map(int, matrix_lines[i].split()))
            for j in range(n):
                if row_vals[j] != 0:
                    k_val = row_vals[j]
                    k_row = i + 1
                    k_col = j + 1
                    found = True
                    break
            if found:
                break
        except (ValueError, IndexError):
            # Handle malformed matrix rows
            continue

    # If k was not found (e.g., all zeros), exit.
    if k_row == -1:
        return

    # The center of an N x N matrix (with 1-based indexing) is at (N/2 + 1, N/2 + 1).
    center_pos = n // 2 + 1

    # The number of moves is the Manhattan distance from the current position to the center.
    # This is because each move (swapping adjacent rows/columns) reduces the distance in one axis by 1.
    row_moves = abs(k_row - center_pos)
    col_moves = abs(k_col - center_pos)
    total_moves = row_moves + col_moves

    # Print the integers k, r, c, and z (total_moves) as requested.
    print(f"{k_val} {k_row} {k_col} {total_moves}")

solve()