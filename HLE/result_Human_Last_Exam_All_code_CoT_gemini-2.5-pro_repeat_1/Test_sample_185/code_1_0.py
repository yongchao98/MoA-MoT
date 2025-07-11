import sys

def solve():
    """
    Solves the beautiful matrix problem.
    """
    try:
        # Read the matrix dimension N
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)
        
        k_val = 0
        k_row = -1
        k_col = -1

        # Read the matrix and find the non-zero digit's position
        for i in range(n):
            row_str = sys.stdin.readline()
            if not row_str:
                return
            # Split the line into numbers
            row_vals = list(map(int, row_str.split()))
            # Check if the non-zero digit is in this row
            if 1 in row_vals or 2 in row_vals or 3 in row_vals or 4 in row_vals or 5 in row_vals or 6 in row_vals or 7 in row_vals or 8 in row_vals or 9 in row_vals:
                for j in range(n):
                    if row_vals[j] != 0:
                        k_val = row_vals[j]
                        # Convert 0-based index to 1-based index
                        k_row = i + 1
                        k_col = j + 1
                        break # Found it in this row
            if k_val != 0:
                # If we've found the number, we don't need to search other rows
                # but we must consume the rest of the input lines.
                pass

        # The center coordinate (1-based)
        center = n // 2 + 1

        # Calculate the Manhattan distance (minimum moves)
        row_moves = abs(k_row - center)
        col_moves = abs(k_col - center)
        total_moves = row_moves + col_moves

        # Print the required output
        # The equation for the result z is: z = abs(r - center) + abs(c - center)
        # The numbers in the equation are r, c, center, and z.
        # The final output requires k, r, c, and z.
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential errors gracefully
        pass

solve()
