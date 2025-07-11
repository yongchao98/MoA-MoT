def solve():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """
    grid_str = """
◫◩▤▤◧◨▥▣▩◨◪
◫◫◫◧◨◪◩▩▩◨▦
▦▧◫▨◧◧◪▥▤⧰◫
◧◫▣◩◫◨▨◪▩▤◨
▦▨◪◪▣▧▩▦◨▩▨
◨◫◫◪◪▨▥◪◧▩◨
◧▨▤◩◫◫▣◫◧◨▥
▩▦▥▩◧◧▧▣◪◨◪
◪◨◫▧◫▩▧◧◩⧰▩
◩▨▤▨◫▩◨◨◨◫▥
▤▨◫◧◨◪▣▤◨▥◧
""".replace("⧰","\u25A7")

    grid = [list(row) for row in grid_str.strip().split('\n')]

    def print_grid(current_grid):
        for row in current_grid:
            print("".join(row))

    # --- Transformation Functions ---
    def rotate_cw(g):
        """Rotates 90 degrees clockwise."""
        return [list(row) for row in zip(*g[::-1])]

    def flip_vertical(g):
        """Flips vertically (upside-down)."""
        return g[::-1]

    def transpose(g):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*g)]

    def rotate_ccw(g):
        """Rotates 90 degrees counter-clockwise."""
        return [list(row) for row in reversed(list(zip(*g)))]

    def rotate_180(g):
        """Rotates 180 degrees."""
        return [row[::-1] for row in g[::-1]]

    def flip_horizontal(g):
        """Flips horizontally (mirror)."""
        return [row[::-1] for row in g]

    # --- Applying the 10-step algorithm ---
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
    
    print_grid(grid)

solve()