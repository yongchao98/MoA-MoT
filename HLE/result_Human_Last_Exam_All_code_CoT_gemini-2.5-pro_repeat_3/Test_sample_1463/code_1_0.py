def solve_grid_puzzle():
    """
    This function solves the grid transformation puzzle by applying a 10-step algorithm.
    """
    # The initial 11x11 grid represented as a list of lists.
    grid_str = """◫◩▤▤◧◨▥▣▩◨◪
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
    
    grid = [list(row) for row in grid_str.split('\n')]

    # --- Transformation Functions ---

    def transpose(g):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def flip_horizontal(g):
        """Flips the grid horizontally (left-to-right mirror)."""
        return [row[::-1] for row in g]

    def rotate_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return flip_horizontal(transpose(g))

    def rotate_counter_clockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        return flip_vertical(transpose(g))

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return flip_horizontal(flip_vertical(g))

    # --- Applying the 10-step Algorithm ---

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

    # --- Print the Final Grid ---
    for row in grid:
        print("".join(row))

solve_grid_puzzle()