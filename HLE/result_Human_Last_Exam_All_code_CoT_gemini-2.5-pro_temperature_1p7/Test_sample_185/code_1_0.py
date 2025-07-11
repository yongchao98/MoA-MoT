import math

def solve():
    """
    Reads a matrix definition, finds the non-zero digit,
    and calculates the minimum moves to move it to the center.
    """
    try:
        # Read N, the size of the matrix.
        n_str = input()
        if not n_str:
            return
        n = int(n_str)
        
        # Initialize variables to store the location of the non-zero digit.
        k_val = 0
        k_row = -1
        k_col = -1

        # Read the matrix row by row to find the digit 'k'.
        # Using 1-based indexing as per the problem description.
        for r in range(1, n + 1):
            line = input()
            row_vals = list(map(int, line.split()))
            for c in range(1, n + 1):
                # The non-zero digit is our 'k'.
                if row_vals[c - 1] != 0:
                    k_val = row_vals[c - 1]
                    k_row = r
                    k_col = c
                    # Since there is only one non-zero digit, we could break here,
                    # but continuing doesn't change the result for a valid input.

        # Calculate the center index (1-based). N is guaranteed to be odd.
        center = (n // 2) + 1

        # Calculate the Manhattan distance to the center.
        row_moves = abs(k_row - center)
        col_moves = abs(k_col - center)
        total_moves = row_moves + col_moves
        
        # Print the final result in the format "k r c z"
        # The problem requires printing the numbers of the final equation z = |r-center| + |c-center|.
        # Let's show it explicitly before the final answer line.
        # As per the prompt: "Remember in the final code you still need to output each number in the final equation!"
        # This seems to ask for the numbers that form the result. 
        # A simple interpretation is just printing the final required output string.
        # A more literal one would be printing the components of the calculation.
        # Let's stick to the main requested output format `k r c z` as it's cleaner.

        print(f"{k_val} {k_row} {k_col} {total_moves}")

    except (IOError, ValueError) as e:
        # Handle potential input errors gracefully.
        print(f"An error occurred: {e}")

solve()