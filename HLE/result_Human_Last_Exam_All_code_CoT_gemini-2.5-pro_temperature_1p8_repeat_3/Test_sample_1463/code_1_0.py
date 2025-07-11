def print_grid(grid):
    """Helper function to print the grid."""
    for row in grid:
        print(row)

def transpose(grid):
    """Transposes the grid (swaps rows and columns)."""
    return tuple("".join(row) for row in zip(*grid))

def flip_v(grid):
    """Flips the grid vertically (upside-down)."""
    return tuple(grid[::-1])

def flip_h(grid):
    """Flips the grid horizontally (left-to-right mirror)."""
    return tuple(row[::-1] for row in grid)

def rotate_90_cw(grid):
    """Rotates the grid 90 degrees clockwise."""
    # This is equivalent to a transpose followed by a horizontal flip.
    return flip_h(transpose(grid))

def rotate_90_ccw(grid):
    """Rotates the grid 90 degrees counterclockwise."""
    # This is equivalent to a transpose followed by a vertical flip.
    return flip_v(transpose(grid))
    
def rotate_180(grid):
    """Rotates the grid 180 degrees."""
    # This is equivalent to a horizontal flip followed by a vertical flip.
    return flip_h(flip_v(grid))

def main():
    """Main function to perform the grid transformations."""
    # Initial 11x11 grid
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

    current_grid = tuple(grid_str.strip().split('\n'))

    # The 10-step algorithm
    # 1. Rotate the grid 90 degrees clockwise.
    current_grid = rotate_90_cw(current_grid)
    # 2. Flip the grid vertically (upside-down).
    current_grid = flip_v(current_grid)
    # 3. Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)
    # 4. Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_90_ccw(current_grid)
    # 5. Flip the grid vertically (upside-down).
    current_grid = flip_v(current_grid)
    # 6. Rotate the grid 180 degrees.
    current_grid = rotate_180(current_grid)
    # 7. Flip the grid horizontally (mirror it left-to-right).
    current_grid = flip_h(current_grid)
    # 8. Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)
    # 9. Flip the grid vertically (upside-down).
    current_grid = flip_v(current_grid)
    # 10. Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_90_ccw(current_grid)

    # Print the final grid
    print_grid(current_grid)

if __name__ == "__main__":
    main()