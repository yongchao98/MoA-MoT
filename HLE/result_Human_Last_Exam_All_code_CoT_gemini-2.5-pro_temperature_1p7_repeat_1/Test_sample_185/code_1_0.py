import sys

def solve():
    """
    Reads a matrix definition, finds the non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # 1. Read the Matrix Size (N)
        n = int(sys.stdin.readline())
        
        # 2. Find the Non-Zero Digit and its location
        k_val = 0
        k_row = -1
        k_col = -1
        
        center_pos = n // 2 + 1
        
        for r in range(1, n + 1):
            line = sys.stdin.readline()
            # If we've already found the number, just read the line to advance the input stream
            if k_row != -1:
                continue
                
            row_vals = list(map(int, line.split()))
            for c in range(1, n + 1):
                # The list is 0-indexed, so we access with c-1
                if row_vals[c-1] != 0:
                    # 3. Store Location and Value
                    k_val = row_vals[c-1]
                    k_row = r
                    k_col = c
                    # Once found, no need to check the rest of this row
                    break
        
        # 5. Calculate Minimum Moves
        row_moves = abs(k_row - center_pos)
        col_moves = abs(k_col - center_pos)
        total_moves = row_moves + col_moves
        
        # 6. Print the Result
        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential empty input or non-integer input
        print(f"Error processing input: {e}", file=sys.stderr)

solve()