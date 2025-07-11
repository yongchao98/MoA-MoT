import sys

def solve():
    """
    Reads a matrix from stdin, finds the single non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read the size of the matrix
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str)

        k_val = 0
        k_row = 0
        k_col = 0
        found = False

        # Loop through the matrix to find the non-zero digit 'k'
        # We use 1-based indexing for rows and columns as per the problem.
        for i in range(1, n + 1):
            line = sys.stdin.readline()
            if not line:
                continue
            
            # To save a little time, only parse the line if we haven't found the digit yet
            if not found:
                row_values = list(map(int, line.split()))
                for j in range(1, n + 1):
                    # In a list, column j (1-based) is at index j-1 (0-based)
                    if row_values[j - 1] != 0:
                        k_val = row_values[j - 1]
                        k_row = i
                        k_col = j
                        found = True
                        break
        
        # Calculate the center coordinate (1-based)
        center = n // 2 + 1
        
        # Calculate the Manhattan distance to the center
        moves = abs(k_row - center) + abs(k_col - center)
        
        # Print the final result as "k r c z"
        print(k_val, k_row, k_col, moves)

    except (IOError, ValueError) as e:
        # Handle potential errors with input
        # print(f"An error occurred: {e}", file=sys.stderr)
        pass

solve()