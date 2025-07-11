def solve():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """

    # Helper function to print the grid in a readable format.
    def print_grid(grid):
        for row in grid:
            print("".join(row))

    # --- Transformation Functions ---

    def transpose(grid):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*grid)]

    def flip_vertical(grid):
        """Flips the grid upside-down."""
        return grid[::-1]

    def flip_horizontal(grid):
        """Mirrors the grid left-to-right."""
        return [row[::-1] for row in grid]

    def rotate_clockwise(grid):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to transposing and then flipping horizontally.
        return flip_horizontal(transpose(grid))

    def rotate_counter_clockwise(grid):
        """Rotates the grid 90 degrees counter-clockwise."""
        # This is equivalent to transposing and then flipping vertically.
        return flip_vertical(transpose(grid))
    
    def rotate_180(grid):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping vertically and then horizontally.
        return flip_horizontal(flip_vertical(grid))

    # Initial 11x11 grid provided in the problem.
    initial_grid_str = """◫◩▤▤◧◨▥▣▩◨◪
◫◫◫◧◨◪◩▩▩◨▦
▦▧◫▨◧◧◪▥▤▧◫
◧◫▣◩◫◨▨◪▩▤◨
▦▨◪◪▣▧▩▦◨▩▨
◨◫◫◪◪▨▥◪◧▩◨
◧▨▤◩◫◫▣◫◧◨▥
▩▦▥▩◧◧▧▣◪◨◪
◪◨◫▧◫▩▧◧◩▧▩
◩▨▤▨◫▩◨◨◨◫▥
▤▨◫◧◨◪▣▤◨▥◧"""

    # Convert the string representation into a list of lists of characters.
    grid = [list(row) for row in initial_grid_str.split('\n')]

    # --- Apply the 10-step algorithm ---

    # 1. Rotate the grid 90 degrees clockwise.
    grid = rotate_clockwise(grid)
    
    # 2. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    
    # 3. Transpose the grid (swap rows and columns).
    grid = transpose(grid)
    
    # 4. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_counter_clockwise(grid)
    
    # 5. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    
    # 6. Rotate the grid 180 degrees.
    grid = rotate_180(grid)
    
    # 7. Flip the grid horizontally (mirror it left-to-right).
    grid = flip_horizontal(grid)
    
    # 8. Transpose the grid (swap rows and columns).
    grid = transpose(grid)
    
    # 9. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    
    # 10. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_counter_clockwise(grid)

    # Print the final grid to the console.
    print_grid(grid)

solve()