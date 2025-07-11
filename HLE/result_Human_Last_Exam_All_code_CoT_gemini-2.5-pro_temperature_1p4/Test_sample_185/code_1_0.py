import sys
import math

def solve():
    """
    Reads a matrix from standard input, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read the size of the matrix
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())
        
        # Calculate the center coordinate (1-based index)
        center = n // 2 + 1
        
        found_k = 0
        found_r = 0
        found_c = 0

        # Loop through rows to find the non-zero number
        for i in range(1, n + 1):
            line = sys.stdin.readline().strip()
            if not line:
                continue
            row_vals = list(map(int, line.split()))
            
            # Check if the non-zero number is in this row
            if any(val != 0 for val in row_vals):
                for j in range(1, n + 1):
                    if row_vals[j-1] != 0:
                        found_k = row_vals[j-1]
                        found_r = i
                        found_c = j
                        # Since there is only one non-zero number, we can stop
                        break
            if found_k != 0:
                break

        # Calculate the minimum number of moves (Manhattan distance)
        moves = abs(found_r - center) + abs(found_c - center)

        # Print the final result: k, original_row, original_col, moves
        print(f"{found_k} {found_r} {found_c} {moves}")

    except (IOError, ValueError) as e:
        # Handle potential input errors
        print(f"An error occurred: {e}", file=sys.stderr)

solve()