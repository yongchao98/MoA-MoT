import sys

def solve():
    """
    Reads a matrix definition from stdin, finds the single non-zero element,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read N, the size of the matrix.
        # The first line is N. Subsequent N lines are the matrix rows.
        # Reading all lines to properly handle input stream.
        lines = sys.stdin.readlines()
        n_str = lines[0]
        n = int(n_str.strip())
        
        matrix_lines = lines[1:n+1]

        k_val = 0
        k_r, k_c = 0, 0

        # Find the non-zero number 'k' and its position (r, c) using 1-based indexing.
        for i in range(n):
            row_str = matrix_lines[i]
            row = list(map(int, row_str.split()))
            for j in range(n):
                if row[j] != 0:
                    k_val = row[j]
                    k_r = i + 1
                    k_c = j + 1
                    # Since there's only one non-zero digit, we could break here,
                    # but iterating through all is fine for the given constraints
                    # and ensures all input is consumed.

        # The center of the matrix (using 1-based indexing)
        center = n // 2 + 1

        # The minimum number of moves is the Manhattan distance to the center.
        moves = abs(k_r - center) + abs(k_c - center)

        # Print the required output: k, original row, original column, and minimum moves.
        print(k_val, k_r, k_c, moves)

    except (IOError, ValueError) as e:
        # Handle potential errors in input reading or conversion
        # This part is for robustness and not strictly required by the problem statement
        pass

solve()