def solve_grid_puzzle():
    """
    This function applies a 10-step transformation algorithm to an 11x11 grid
    and prints the final result.
    """

    def rotate_cw(grid):
        """Rotates the grid 90 degrees clockwise."""
        return [list(row[::-1]) for row in zip(*grid)]

    def flip_vertical(grid):
        """Flips the grid vertically (upside-down)."""
        return grid[::-1]

    def transpose(grid):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*grid)]

    def rotate_ccw(grid):
        """Rotates the grid 90 degrees counter-clockwise."""
        transposed = [list(row) for row in zip(*grid)]
        return transposed[::-1]

    def rotate_180(grid):
        """Rotates the grid 180 degrees."""
        return [row[::-1] for row in grid[::-1]]

    def flip_horizontal(grid):
        """Flips the grid horizontally (mirror image)."""
        return [row[::-1] for row in grid]

    # The initial 11x11 grid.
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

    grid = [list(row) for row in initial_grid_str.strip().split('\n')]

    # Apply the 10-step algorithm sequentially.
    # 1. Rotate the grid 90 degrees clockwise.
    grid = rotate_cw(grid)
    # 2. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)
    # 3. Transpose the grid (swap rows and columns).
    grid = transpose(grid)
    # 4. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_ccw(grid)
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
    grid = rotate_ccw(grid)

    # Print the final state of the grid.
    for row in grid:
        print("".join(row))

solve_grid_puzzle()