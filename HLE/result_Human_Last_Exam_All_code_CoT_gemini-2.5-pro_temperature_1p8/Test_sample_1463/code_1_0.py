import copy

def main():
    """
    Solves the grid transformation puzzle by applying a 10-step algorithm.
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

    grid = [list(row) for row in grid_str.split('\n')]

    def print_grid(g, title):
        """Helper function to print the grid with a title."""
        # This function is for debugging and not used in the final output
        # to adhere to the single code block rule.
        pass

    def transpose(g):
        """Swaps rows and columns of the grid."""
        return [list(row) for row in zip(*g)]

    def flip_vertical(g):
        """Flips the grid upside-down."""
        return g[::-1]

    def flip_horizontal(g):
        """Mirrors the grid left-to-right."""
        return [row[::-1] for row in g]

    def rotate_clockwise(g):
        """Rotates the grid 90 degrees clockwise."""
        return flip_horizontal(transpose(g))

    def rotate_counter_clockwise(g):
        """Rotates the grid 90 degrees counterclockwise."""
        return flip_vertical(transpose(g))
        
    def rotate_180(g):
        """Rotates the grid 180 degrees."""
        return flip_horizontal(flip_vertical(g))

    # Apply the 10-step algorithm
    current_grid = copy.deepcopy(grid)
    
    # 1. Rotate the grid 90 degrees clockwise.
    current_grid = rotate_clockwise(current_grid)
    
    # 2. Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)
    
    # 3. Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)
    
    # 4. Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_counter_clockwise(current_grid)
    
    # 5. Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)
    
    # 6. Rotate the grid 180 degrees.
    current_grid = rotate_180(current_grid)
    
    # 7. Flip the grid horizontally (mirror it left-to-right).
    current_grid = flip_horizontal(current_grid)
    
    # 8. Transpose the grid (swap rows and columns).
    current_grid = transpose(current_grid)
    
    # 9. Flip the grid vertically (upside-down).
    current_grid = flip_vertical(current_grid)
    
    # 10. Rotate the grid 90 degrees counterclockwise.
    current_grid = rotate_counter_clockwise(current_grid)

    # Print the final grid
    for row in current_grid:
        print("".join(row))

if __name__ == "__main__":
    main()