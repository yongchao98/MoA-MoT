def solve_grid_puzzle():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
    """

    # Helper function to print the grid
    def print_grid(grid, step_name):
        # This is for debugging and showing intermediate steps if needed.
        # For the final answer, we only need to print the last grid.
        pass

    # --- Transformation Functions ---

    def rotate_90_clockwise(grid):
        """Rotates the grid 90 degrees clockwise."""
        # This is equivalent to transposing and then reversing each row.
        return [list(row[::-1]) for row in zip(*grid)]

    def flip_vertical(grid):
        """Flips the grid vertically (upside-down)."""
        return grid[::-1]

    def transpose(grid):
        """Transposes the grid (swaps rows and columns)."""
        return [list(row) for row in zip(*grid)]

    def rotate_90_counterclockwise(grid):
        """Rotates the grid 90 degrees counterclockwise."""
        # This is equivalent to transposing and then reversing the order of the new rows.
        transposed = [list(row) for row in zip(*grid)]
        return transposed[::-1]

    def rotate_180(grid):
        """Rotates the grid 180 degrees."""
        # This is equivalent to flipping vertically and then horizontally.
        return [row[::-1] for row in grid[::-1]]

    def flip_horizontal(grid):
        """Flips the grid horizontally (left-to-right mirror)."""
        return [row[::-1] for row in grid]

    # --- Initial Grid ---

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

    # --- Apply 10-Step Algorithm ---
    
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

    # --- Print Final Result ---
    
    for row in grid:
        print("".join(row))

solve_grid_puzzle()