def solve_grid_puzzle():
    """
    Solves the 10-step grid transformation puzzle.
    """
    # 1. Represent the initial 11x11 grid.
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

    grid = [list(row) for row in grid_str.split('\n')]

    # 2. Define the transformation functions.
    def rotate_90_clockwise(g):
        return [list(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        return g[::-1]

    def transpose(g):
        return [list(row) for row in zip(*g)]

    def rotate_90_counterclockwise(g):
        return [list(row) for row in zip(*g)][::-1]

    def rotate_180(g):
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        return [row[::-1] for row in g]

    # 3. Apply the 10-step algorithm sequentially.
    current_grid = grid

    # Step 1: Rotate the grid 90 degrees clockwise.
    current_grid = rotate_90_clockwise(current_grid)

    # Step 2: Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)

    # Step 3: Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)

    # Step 4: Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_90_counterclockwise(current_grid)

    # Step 5: Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)

    # Step 6: Rotate the grid 180 degrees.
    current_grid = rotate_180(current_grid)

    # Step 7: Flip the grid horizontally (mirror it left-to-right).
    current_grid = flip_horizontal(current_grid)

    # Step 8: Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)

    # Step 9: Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)

    # Step 10: Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_90_counterclockwise(current_grid)
    
    # 4. Print the final grid state.
    for row in current_grid:
        print("".join(row))

solve_grid_puzzle()