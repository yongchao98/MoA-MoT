def solve():
    """
    Solves the grid transformation puzzle by applying a series of operations.
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

    # Convert the string grid into a list of lists of characters
    grid = [list(line) for line in grid_str.strip().split('\n')]

    # --- Transformation Functions ---

    def transpose(g):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        """Flips the grid upside-down."""
        return g[::-1]

    def flip_horizontal(g):
        """Flips the grid left-to-right."""
        return [row[::-1] for row in g]

    def rotate_90_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to a transpose followed by a horizontal flip.
        return [list(r)[::-1] for r in zip(*g)]

    def rotate_90_counterclockwise(g):
        """Rotates the grid 90 degrees counterclockwise."""
        # This is equivalent to a transpose followed by a vertical flip.
        return [list(r) for r in zip(*g)][::-1]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        # This is equivalent to a horizontal flip followed by a vertical flip.
        return [row[::-1] for row in g[::-1]]

    # --- Apply the 10-step Algorithm ---

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

    # --- Print the Final Grid ---
    for row in grid:
        print("".join(row))

solve()