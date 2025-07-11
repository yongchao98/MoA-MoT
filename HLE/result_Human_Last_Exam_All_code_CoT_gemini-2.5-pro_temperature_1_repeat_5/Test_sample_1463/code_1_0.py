def print_grid(grid):
    """Prints the grid in a readable format."""
    for row in grid:
        print("".join(row))

def rotate_clockwise(grid):
    """Rotates the grid 90 degrees clockwise."""
    # This is equivalent to transposing and then flipping horizontally.
    return [list(row) for row in zip(*grid[::-1])]

def flip_vertical(grid):
    """Flips the grid vertically (upside-down)."""
    return grid[::-1]

def transpose(grid):
    """Transposes the grid (swaps rows and columns)."""
    return [list(row) for row in zip(*grid)]

def rotate_counter_clockwise(grid):
    """Rotates the grid 90 degrees counter-clockwise."""
    # This is equivalent to transposing and then flipping vertically.
    return [list(row) for row in reversed(list(zip(*grid)))]

def rotate_180(grid):
    """Rotates the grid 180 degrees."""
    # This is equivalent to flipping vertically and then horizontally.
    return [row[::-1] for row in grid[::-1]]

def flip_horizontal(grid):
    """Flips the grid horizontally (mirror image)."""
    return [row[::-1] for row in grid]

def main():
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

    grid = [list(row) for row in grid_str.split('\n')]

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

    # Print the final grid
    print_grid(grid)

if __name__ == "__main__":
    main()