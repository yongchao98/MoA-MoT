def solve_grid_puzzle():
    """
    This function applies a 10-step algorithm to a grid and prints the final state.
    """
    # Initial 11x11 grid
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

    def rotate_90_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return [''.join(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [''.join(col) for col in zip(*g)]

    def rotate_90_counterclockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        transposed = list(zip(*g))
        return [''.join(col) for col in transposed[::-1]]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        """Flips the grid horizontally (left-to-right mirror)."""
        return [row[::-1] for row in g]

    # --- Applying the 10-step algorithm ---

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
        print(row)

solve_grid_puzzle()