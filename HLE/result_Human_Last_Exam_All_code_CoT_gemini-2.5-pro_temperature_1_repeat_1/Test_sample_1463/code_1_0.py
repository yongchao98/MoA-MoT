import sys

def solve():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """

    # The initial 11x11 grid as a list of lists of characters.
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
    grid = [list(row) for row in initial_grid_str.strip().split('\n')]

    # --- Transformation Functions ---

    def rotate_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return [list(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips the grid vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*g)]

    def rotate_counter_clockwise(g):
        """Rotates the grid 90 degrees counterclockwise."""
        # This is equivalent to transposing and then flipping vertically.
        return [list(row) for row in zip(*g)][::-1]

    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping vertically and then horizontally.
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        """Flips the grid horizontally (mirror image)."""
        return [row[::-1] for row in g]

    # --- Apply the 10-step Algorithm ---

    # 1. Rotate the grid 90 degrees clockwise.
    grid = rotate_clockwise(grid)

    # 2. Flip the grid vertically (upside-down).
    grid = flip_vertical(grid)

    # 3. Transpose the grid (swap rows and columns).
    grid = transpose(grid)

    # 4. Rotate the grid 90 degrees counterclockwise.
    grid = rotate_counter_clockwise(grid)

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
    grid = rotate_counter_clockwise(grid)

    # --- Print the Final Result ---
    for row in grid:
        print("".join(row))

solve()