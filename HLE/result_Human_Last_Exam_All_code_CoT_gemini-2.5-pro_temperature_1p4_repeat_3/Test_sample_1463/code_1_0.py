def solve_grid_puzzle():
    """
    This function applies a 10-step transformation algorithm to an 11x11 grid
    and prints the final result.
    """
    # The initial 11x11 grid represented as a list of strings.
    grid = [
        "◫◩▤▤◧◨▥▣▩◨◪",
        "◫◫◫◧◨◪◩▩▩◨▦",
        "▦▧◫▨◧◧◪▥▤▧◫",
        "◧◫▣◩◫◨▨◪▩▤◨",
        "▦▨◪◪▣▧▩▦◨▩▨",
        "◨◫◫◪◪▨▥◪◧▩◨",
        "◧▨▤◩◫◫▣◫◧◨▥",
        "▩▦▥▩◧◧▧▣◪◨◪",
        "◪◨◫▧◫▩▧◧◩▧▩",
        "◩▨▤▨◫▩◨◨◨◫▥",
        "▤▨◫◧◨◪▣▤◨▥◧",
    ]

    # --- Transformation Functions ---

    def rotate_90_cw(g):
        """Rotates the grid 90 degrees clockwise."""
        return [''.join(c) for c in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [''.join(c) for c in zip(*g)]

    def rotate_90_ccw(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        transposed_grid = list(zip(*g))
        return [''.join(row) for row in transposed_grid[::-1]]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        """Flips the grid horizontally (mirror image)."""
        return [row[::-1] for row in g]

    # --- Applying the 10-step algorithm ---

    # 1. Rotate the grid 90 degrees clockwise.
    grid = rotate_90_cw(grid)
    
    # 2. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    
    # 3. Transpose the grid (swap rows and columns).
    grid = transpose(grid)
    
    # 4. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_90_ccw(grid)
    
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
    grid = rotate_90_ccw(grid)

    # Print the final state of the grid.
    for row in grid:
        print(row)

solve_grid_puzzle()