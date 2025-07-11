def rotate_90_clockwise(grid):
    """Rotates the grid 90 degrees clockwise."""
    return [list(reversed(col)) for col in zip(*grid)]

def rotate_90_counterclockwise(grid):
    """Rotates the grid 90 degrees counterclockwise."""
    return [list(col) for col in reversed(list(zip(*grid)))]

def rotate_180(grid):
    """Rotates the grid 180 degrees."""
    return [row[::-1] for row in grid[::-1]]

def flip_vertical(grid):
    """Flips the grid vertically (upside-down)."""
    return grid[::-1]

def flip_horizontal(grid):
    """Flips the grid horizontally (mirror image)."""
    return [row[::-1] for row in grid]

def transpose(grid):
    """Transposes the grid (swaps rows and columns)."""
    return [list(col) for col in zip(*grid)]

def solve():
    """
    Applies the 10-step algorithm to the grid and prints the final result.
    """
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

    # The 10-step algorithm
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
    for row in grid:
        print("".join(row))

solve()