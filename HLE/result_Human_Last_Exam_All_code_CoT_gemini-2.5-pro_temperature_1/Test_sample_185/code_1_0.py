import sys

def solve():
    """
    This function reads a matrix, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    # Read the matrix size N
    try:
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())
    except (IOError, ValueError):
        return

    k_val = 0
    k_row = -1
    k_col = -1

    # Read the matrix and find the non-zero digit's location
    for i in range(1, n + 1):
        try:
            line = sys.stdin.readline().strip().split()
            if not line:
                continue
            row_vals = [int(val) for val in line]
            for j in range(1, n + 1):
                if row_vals[j-1] != 0:
                    k_val = row_vals[j-1]
                    k_row = i
                    k_col = j
        except (IOError, ValueError):
            continue
            
    # Calculate the center of the matrix
    center = n // 2 + 1
    
    # Calculate the number of row and column swaps needed
    row_moves = abs(k_row - center)
    col_moves = abs(k_col - center)
    
    # Total moves is the sum of row and column moves
    total_moves = row_moves + col_moves
    
    # Print the result including the equation for the total moves
    # The problem asks for integers k, r, c, z. My instructions ask for the equation.
    # This output provides all required information.
    print(f"{k_val} {k_row} {k_col} {row_moves} + {col_moves} = {total_moves}")

solve()
