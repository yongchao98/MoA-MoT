import sys

def solve():
    """
    Reads a matrix, finds a single non-zero digit, and calculates
    the minimum moves to move it to the center.
    """
    try:
        # Read the size of the matrix
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle case of empty input at the end of a file
        return

    # Initialize variables to store the location of the non-zero digit
    k_value = 0
    k_row = 0
    k_col = 0

    # Read the matrix and find the coordinates of the non-zero digit
    for i in range(1, n + 1):
        # Read a row and split it into numbers
        row_values = list(map(int, sys.stdin.readline().strip().split()))
        
        # Check if the non-zero digit is in this row
        if 1 in row_values or 2 in row_values or 3 in row_values or 4 in row_values or 5 in row_values or 6 in row_values or 7 in row_values or 8 in row_values or 9 in row_values:
            for j in range(len(row_values)):
                if row_values[j] != 0:
                    k_value = row_values[j]
                    k_row = i
                    k_col = j + 1 # Convert 0-based index to 1-based
                    break # Found it, no need to check rest of the row

    # The center coordinate for an N x N matrix (where N is odd)
    center_coord = n // 2 + 1
    
    # Calculate the number of moves (Manhattan distance)
    # Moves = |current_row - center_row| + |current_col - center_col|
    moves_needed = abs(k_row - center_coord) + abs(k_col - center_coord)
    
    # Print the result in the specified format: k r c z
    print(f"{k_value} {k_row} {k_col} {moves_needed}")

solve()