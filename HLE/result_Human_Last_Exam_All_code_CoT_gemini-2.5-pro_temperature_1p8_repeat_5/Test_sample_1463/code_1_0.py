def solve_grid_puzzle():
    """
    Applies a 10-step transformation algorithm to an 11x11 grid and prints the result.
    """

    # The initial 11x11 grid
    initial_grid_str = [
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

    # Convert the grid of strings to a list of lists of characters for easier manipulation
    grid = [list(row) for row in initial_grid_str]

    # --- Transformation Functions ---

    def rotate_90_clockwise(g):
        # Transpose and then flip horizontally
        return [list(row)[::-1] for row in zip(*g)]

    def flip_vertical(g):
        # Reverse the order of rows
        return g[::-1]

    def transpose(g):
        # Swap rows and columns
        return [list(col) for col in zip(*g)]

    def rotate_90_counterclockwise(g):
        # Transpose and then flip vertically
        transposed = list(zip(*g))
        return [list(row) for row in transposed[::-1]]

    def rotate_180(g):
        # Flip vertically and then horizontally
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        # Reverse the order of characters in each row
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

    # Print the final grid configuration
    for row in grid:
        print("".join(row))

solve_grid_puzzle()