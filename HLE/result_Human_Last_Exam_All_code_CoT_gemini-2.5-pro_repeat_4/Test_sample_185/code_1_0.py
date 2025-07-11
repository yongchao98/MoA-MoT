import sys
import math

def solve_beautiful_matrix():
    """
    Reads a matrix from stdin, finds the single non-zero digit, and calculates
    the minimum number of swaps to move it to the center.
    """
    # Initialize variables to store the location and value of the non-zero digit.
    k_val = 0
    k_row = -1
    k_col = -1
    
    # Read all lines from standard input.
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
            
        # The first line contains N.
        n = int(lines[0])
        
        # The center coordinate (using 1-based indexing).
        center = n // 2 + 1
        
        # Iterate through the matrix data (lines 1 to N+1) to find the digit.
        for i in range(n):
            # Read a row of numbers.
            row_data = list(map(int, lines[i + 1].split()))
            for j in range(n):
                # Check if the current element is the non-zero digit.
                if row_data[j] != 0:
                    k_val = row_data[j]
                    # Store its 1-based row and column index.
                    k_row = i + 1
                    k_col = j + 1
                    # Since there is only one such digit, we can stop searching.
                    break
            if k_val != 0:
                break
                
        # Calculate the minimum number of moves (Manhattan distance).
        moves = abs(k_row - center) + abs(k_col - center)
        
        # Print the final result as required: k r c z
        print(f"{k_val} {k_row} {k_col} {moves}")

    except (IOError, ValueError, IndexError):
        # Handle cases with invalid input format.
        # For this problem, we assume the input is always valid.
        pass

solve_beautiful_matrix()