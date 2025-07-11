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
    
    grid = [list(row) for row in grid_str.splitlines()]

    # --- Transformation Functions ---
    
    def rotate_90_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to transposing and then reversing each new row.
        # Or, zip of the reversed original rows.
        return [list(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*g)]

    def rotate_90_counterclockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        # This is equivalent to transposing and then reversing the new list of rows.
        return [list(row) for row in reversed(list(zip(*g)))]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        # Flip vertically, then horizontally.
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

    # --- Print the final grid ---
    for row in grid:
        print("".join(row))

solve()