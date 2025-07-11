def solve_grid_puzzle():
    """
    This function applies a 10-step algorithm to a grid and prints the final state.
    """
    # The initial 11x11 grid represented as a list of lists of characters.
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
        list("▤▨◫◧◨◪▣▤◨▥◧")
    ]

    # --- Transformation Functions ---

    def rotate_90_clockwise(g):
        # Transpose by zipping and then reverse each new row
        return [list(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        # Reverse the order of rows
        return g[::-1]

    def transpose(g):
        # Swap rows and columns using zip
        return [list(col) for col in zip(*g)]

    def rotate_90_counterclockwise(g):
        # Transpose and then reverse the order of the new rows
        transposed = [list(col) for col in zip(*g)]
        return transposed[::-1]

    def rotate_180(g):
        # Reverse the order of rows, then reverse each individual row
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        # Reverse each individual row
        return [row[::-1] for row in g]

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

    # Print the final grid state
    for row in grid:
        print("".join(row))

solve_grid_puzzle()