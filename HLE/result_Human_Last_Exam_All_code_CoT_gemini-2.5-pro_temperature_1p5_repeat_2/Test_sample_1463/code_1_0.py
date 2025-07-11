def solve():
    """
    Solves the grid transformation puzzle.
    """
    initial_grid_str = """
◫◩▤▤◧◨▥▣▩◨◪
◫◫◫◧◨◪◩▩▩◨▦
▦▧◫▨◧◧◪▥▤▧◫
◧◫▣◩◫◨▨◪▩▤◨
▦▨◪◪▣▧▩▦◨▩▨
◨◫◫◪◪▨▥◪◧▩◨
◧▨▤◩◫◫▣◫◧◨▥
▩▦▥▩◧◧▧▣◪◨◪
◪◨◫▧◫▩▧◧◩▧▩
◩▨▤▨◫▩◨◨◨◫▥
▤▨◫◧◨◪▣▤◨▥◧
"""
    # Parse the string into a list of lists
    grid = [list(row) for row in initial_grid_str.strip().split('\n')]

    # --- Transformation Functions ---
    def rotate_90_clockwise(g):
        # Transpose and then reverse each row
        return [list(row[::-1]) for row in zip(*g)]

    def flip_vertical(g):
        # Reverse the order of rows
        return g[::-1]

    def transpose(g):
        # Swap rows and columns
        return [list(row) for row in zip(*g)]

    def rotate_90_counterclockwise(g):
        # Transpose and then reverse the order of the new rows
        transposed = [list(row) for row in zip(*g)]
        return transposed[::-1]
    
    def rotate_180(g):
        # Reverse each row, then reverse the order of rows
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        # Reverse each row
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

    # --- Print the final result ---
    for row in grid:
        print("".join(row))

solve()