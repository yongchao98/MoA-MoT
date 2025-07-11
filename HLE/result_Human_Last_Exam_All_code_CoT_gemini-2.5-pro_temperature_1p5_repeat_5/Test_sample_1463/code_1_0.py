def solve_grid_puzzle():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """

    # The initial 11x11 grid as a list of lists of characters
    grid = [
        list("◫◩▤▤◧◨▥▣▩◨◪"),
        list("◫◫◫◧◨◪◩▩▩◨▦"),
        list("▦▧◫▨◧◧◪▥▤▧◫"),
        list("◧◫▣◩◫◨▨◪▩▤◨"),
        list("▦▨◪◪▣▧▩▦◨▩▨"),
        list("◨◫◫◪◪▨▥◪◧▩◨"),
        list("◧▨▤◩◫◫▣◫◧◨▥"),
        list("▩▦▥▩◧◧▧▣◪◨◪"),
        list("◪◨◫▧◫▩▧◧◩▧▩"),
        list("◩▨▤▨◫▩◨◨◨◫▥"),
        list("▤▨◫◧◨◪▣▤◨▥◧"),
    ]

    # --- Transformation Functions ---

    def transpose(current_grid):
        """Swaps rows and columns."""
        return [list(row) for row in zip(*current_grid)]

    def flip_vertical(current_grid):
        """Flips the grid upside-down."""
        return current_grid[::-1]

    def flip_horizontal(current_grid):
        """Mirrors the grid left-to-right."""
        return [row[::-1] for row in current_grid]

    def rotate_90_clockwise(current_grid):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to transposing and then flipping horizontally.
        return flip_horizontal(transpose(current_grid))

    def rotate_90_counterclockwise(current_grid):
        """Rotates the grid 90 degrees counterclockwise."""
        # This is equivalent to transposing and then flipping vertically.
        return flip_vertical(transpose(current_grid))

    def rotate_180(current_grid):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping vertically and then horizontally.
        return flip_horizontal(flip_vertical(current_grid))

    # --- Applying the 10-step Algorithm ---

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

    # Print the final grid
    for row in grid:
        print("".join(row))

solve_grid_puzzle()