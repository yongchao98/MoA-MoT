def solve_grid_puzzle():
    """
    This function initializes the grid, applies a 10-step transformation algorithm,
    and prints the final result.
    """
    # The initial 11x11 grid
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

    grid = [list(line) for line in initial_grid_str.strip().split('\n')]

    # --- Transformation Functions ---

    def rotate_90_cw(g):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to transposing and then flipping horizontally.
        return [list(row)[::-1] for row in zip(*g)]

    def flip_v(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose_grid(g):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*g)]

    def rotate_90_ccw(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        # This is equivalent to transposing and then flipping vertically.
        transposed = [list(row) for row in zip(*g)]
        return transposed[::-1]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping vertically and then horizontally.
        return [row[::-1] for row in g[::-1]]

    def flip_h(g):
        """Flips the grid horizontally (left-to-right mirror)."""
        return [row[::-1] for row in g]

    # --- Applying the 10-Step Algorithm ---

    # Step 1: Rotate the grid 90 degrees clockwise.
    grid = rotate_90_cw(grid)

    # Step 2: Flip the grid vertically (upside-down).
    grid = flip_v(grid)

    # Step 3: Transpose the grid (swap rows and columns).
    grid = transpose_grid(grid)

    # Step 4: Rotate the grid 90 degrees counterclockwise.
    grid = rotate_90_ccw(grid)

    # Step 5: Flip the grid vertically (upside-down).
    grid = flip_v(grid)

    # Step 6: Rotate the grid 180 degrees.
    grid = rotate_180(grid)

    # Step 7: Flip the grid horizontally (mirror it left-to-right).
    grid = flip_h(grid)

    # Step 8: Transpose the grid (swap rows and columns).
    grid = transpose_grid(grid)

    # Step 9: Flip the grid vertically (upside-down).
    grid = flip_v(grid)

    # Step 10: Rotate the grid 90 degrees counterclockwise.
    grid = rotate_90_ccw(grid)

    # --- Print the Final Grid ---
    for row in grid:
        # Print each character in the final grid configuration
        print("".join(row))

# Execute the function to see the result
solve_grid_puzzle()