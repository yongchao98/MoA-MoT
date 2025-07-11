import sys

def solve():
    """
    This function solves the grid transformation puzzle by applying a 10-step algorithm.
    """
    # The initial 11x11 grid provided in the problem.
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

    # Convert the string grid into a list of lists of characters.
    grid = [list(row) for row in initial_grid_str.split('\n')]

    # Helper functions for grid transformations.
    def transpose(g):
        # Swap rows and columns.
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        # Reverse the order of rows.
        return g[::-1]

    def flip_horizontal(g):
        # Reverse each row.
        return [row[::-1] for row in g]

    def rotate_clockwise(g):
        # A 90-degree clockwise rotation is equivalent to transposing and then flipping horizontally.
        return flip_horizontal(transpose(g))

    def rotate_counter_clockwise(g):
        # A 90-degree counter-clockwise rotation is equivalent to transposing and then flipping vertically.
        return flip_vertical(transpose(g))

    def rotate_180(g):
        # A 180-degree rotation is equivalent to flipping both vertically and horizontally.
        return flip_horizontal(flip_vertical(g))
        
    # The 10-step algorithm.
    transformations = [
        rotate_clockwise,          # 1. Rotate 90 degrees clockwise
        flip_vertical,             # 2. Flip vertically
        transpose,                 # 3. Transpose
        rotate_counter_clockwise,  # 4. Rotate 90 degrees counterclockwise
        flip_vertical,             # 5. Flip vertically
        rotate_180,                # 6. Rotate 180 degrees
        flip_horizontal,           # 7. Flip horizontally
        transpose,                 # 8. Transpose
        flip_vertical,             # 9. Flip vertically
        rotate_counter_clockwise   # 10. Rotate 90 degrees counterclockwise
    ]

    # Apply each transformation sequentially.
    for transform in transformations:
        grid = transform(grid)

    # Print the final grid configuration.
    for row in grid:
        print("".join(row))

solve()