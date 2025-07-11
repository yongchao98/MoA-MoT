def solve():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
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
    # Parse the initial grid string into a list of lists of characters
    grid = [list(row) for row in initial_grid_str.strip().split('\n')]

    def print_grid(g):
        for row in g:
            print("".join(row))

    def transpose(g):
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        return g[::-1]

    def flip_horizontal(g):
        return [row[::-1] for row in g]

    def rotate_90_clockwise(g):
        return flip_horizontal(transpose(g))

    def rotate_90_counterclockwise(g):
        return flip_vertical(transpose(g))

    def rotate_180(g):
        return flip_horizontal(flip_vertical(g))

    # Apply the 10-step algorithm
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
    print_grid(grid)

solve()