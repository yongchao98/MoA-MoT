def solve():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """
    # The initial 11x11 grid
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

    # Parse the string into a list of lists of characters
    grid = [list(row) for row in grid_str.strip().split('\n')]

    # --- Transformation Functions ---

    def rotate_90_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return [list(t) for t in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [list(t) for t in zip(*g)]

    def rotate_90_counterclockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        # This is equivalent to transposing and then reversing rows.
        return [list(r) for r in reversed(list(zip(*g)))]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        """Flips the grid horizontally (mirror image)."""
        return [row[::-1] for row in g]

    # --- Apply the 10-step algorithm ---

    # 1. Rotate the grid 90 degrees clockwise.
    grid = rotate_90_clockwise(grid)
    
    # 2. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    
    # 3. Transpose the grid (swap rows and columns).
    grid = transpose(grid)
    
    # 4. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_90_counterclockwise(grid)
    
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
    grid = rotate_90_counterclockwise(grid)

    # Print the final state of the grid
    for row in grid:
        print("".join(row))

solve()