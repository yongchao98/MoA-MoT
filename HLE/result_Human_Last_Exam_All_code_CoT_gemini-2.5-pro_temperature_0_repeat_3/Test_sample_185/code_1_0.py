import sys

def solve_matrix_problem():
    """
    Reads a matrix from stdin, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read the matrix size N
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())

        k_val = 0
        k_row = -1
        k_col = -1

        # Find the location of the non-zero digit
        for i in range(n):
            line = sys.stdin.readline()
            if not line:
                return
            # Read the current row and split it into numbers
            row_values = list(map(int, line.strip().split()))
            for j in range(n):
                if row_values[j] != 0:
                    k_val = row_values[j]
                    # Convert 0-based index to 1-based index for output
                    k_row = i + 1
                    k_col = j + 1
                    # Found the number, no need to search further in this row
                    break
            # If we found the number, we could break here too, but reading all input is safer
            # if k_row != -1:
            #     break

        # The center of an N x N matrix (where N is odd) is at (N/2 + 1, N/2 + 1)
        # Using integer division //
        center = n // 2 + 1

        # The number of moves is the Manhattan distance from the current position to the center
        moves = abs(k_row - center) + abs(k_col - center)

        # Print the required output
        print(f"{k_val} {k_row} {k_col} {moves}")

    except (ValueError, IndexError):
        # Handle potential empty lines or malformed input
        pass

solve_matrix_problem()