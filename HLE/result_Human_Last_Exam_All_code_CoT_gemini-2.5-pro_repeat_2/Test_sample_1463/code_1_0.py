def solve():
    """
    Solves the grid transformation puzzle by simulating the 10-step algorithm.
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
    
    # Convert the string grid to a list of lists
    grid = [list(line) for line in grid_str.split('\n')]

    # --- Transformation Functions ---

    def transpose(g):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def flip_horizontal(g):
        """Flips the grid horizontally (mirror image)."""
        return [row[::-1] for row in g]

    def rotate_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to flipping vertically, then transposing.
        return transpose(flip_vertical(g))

    def rotate_counter_clockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        # This is equivalent to flipping horizontally, then transposing.
        return transpose(flip_horizontal(g))

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping horizontally, then vertically.
        return flip_vertical(flip_horizontal(g))

    # --- The 10-Step Algorithm ---

    # A list of the operations to be performed in order
    operations = [
        rotate_clockwise,          # 1. Rotate the grid 90 degrees clockwise.
        flip_vertical,             # 2. Flip the grid vertically (upside-down).
        transpose,                 # 3. Transpose the grid (swap rows and columns).
        rotate_counter_clockwise,  # 4. Rotate the grid 90 degrees counterclockwise.
        flip_vertical,             # 5. Flip the grid vertically (upside-down).
        rotate_180,                # 6. Rotate the grid 180 degrees.
        flip_horizontal,           # 7. Flip the grid horizontally (mirror it left-to-right).
        transpose,                 # 8. Transpose the grid (swap rows and columns).
        flip_vertical,             # 9. Flip the grid vertically (upside-down).
        rotate_counter_clockwise   # 10. Rotate the grid 90 degrees counterclockwise.
    ]

    # Apply each operation sequentially to the grid
    current_grid = grid
    for op in operations:
        current_grid = op(current_grid)

    # Print the final grid configuration
    for row in current_grid:
        print("".join(row))

solve()