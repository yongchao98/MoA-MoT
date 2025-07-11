def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any 5-in-a-row configurations.
    """
    rows = 7
    cols = 8
    
    # We use a coloring scheme f(r, c) = (2*(r-1) + (c-1)) mod 5.
    # Any line of 5 will have one square of each color {0,1,2,3,4}.
    # By removing all squares of one color, we break all lines of 5.
    # We find the color that appears the least number of times to minimize removals.
    
    color_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    
    # Iterate through each square on the board (using 1-based indexing)
    for r in range(1, rows + 1):
        for c in range(1, cols + 1):
            # Calculate the color of the square
            color = (2 * (r - 1) + (c - 1)) % 5
            color_counts[color] += 1
            
    # The minimum number of pieces to remove is the size of the smallest color class.
    min_removals = min(color_counts.values())
    
    print(f"Board dimensions: {rows}x{cols}")
    print("Coloring function: f(r, c) = (2*(r-1) + (c-1)) mod 5")
    print("Number of squares for each color:")
    for color, count in color_counts.items():
        print(f"  Color {color}: {count} squares")
        
    print(f"\nThe minimum number of pieces that must be removed is the smallest count.")
    print(f"Result: {min_removals}")

solve_chess_problem()