import sys

def solve():
    """
    Reads a matrix description, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read the size of the matrix, N.
        n_line = sys.stdin.readline()
        if not n_line:
            return
        n = int(n_line)
    except (ValueError, IndexError):
        # This handles cases with malformed input for N, though the problem guarantees valid input.
        return

    # Initialize variables to store the location and value of the non-zero digit.
    k_val = 0
    r_pos = -1
    c_pos = -1

    # Loop through the matrix to find the non-zero element.
    # The problem uses 1-based indexing, so we'll use 'i' from 1 to N.
    for i in range(1, n + 1):
        # Read a row of the matrix and convert its elements to integers.
        row_vals = list(map(int, sys.stdin.readline().split()))
        
        # Check each element in the row.
        # The column index 'j' is 0-based, so the actual column 'c' is j+1.
        for j in range(n):
            if row_vals[j] != 0:
                k_val = row_vals[j]
                r_pos = i
                c_pos = j + 1
                # Since there is only one non-zero number, we don't need to continue
                # searching after this point. However, we must consume the remaining input.
    
    # The center of an N x N matrix (with 1-based indexing) is at (N/2 + 1, N/2 + 1).
    # Using integer division `//` handles this correctly for odd N.
    target_pos = n // 2 + 1
    
    # The minimum number of moves is the Manhattan distance.
    # It's the sum of the absolute differences of the coordinates.
    moves_needed = abs(r_pos - target_pos) + abs(c_pos - target_pos)
    
    # Print the required output: k, r, c, and z (the number of moves).
    print(k_val, r_pos, c_pos, moves_needed)

solve()