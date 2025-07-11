import sys

def solve():
    """
    This script applies a 10-step transformation algorithm to an 11x11 grid
    and prints the final result.
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
    grid = [list(row) for row in grid_str.strip().split('\n')]

    # --- Transformation Functions ---
    def transpose(g):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        """Flips the grid upside-down."""
        return g[::-1]

    def flip_horizontal(g):
        """Flips the grid left-to-right."""
        return [row[::-1] for row in g]

    def rotate_90_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return flip_horizontal(transpose(g))

    def rotate_90_counterclockwise(g):
        """Rotates the grid 90 degrees counter-clockwise."""
        return flip_vertical(transpose(g))

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return flip_horizontal(flip_vertical(g))

    # --- The 10-step Algorithm ---
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

    # --- Print the final grid ---
    for row in grid:
        print("".join(row))

solve()