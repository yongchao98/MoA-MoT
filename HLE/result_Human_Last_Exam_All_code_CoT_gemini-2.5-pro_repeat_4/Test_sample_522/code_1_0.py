def solve_puzzle():
    """
    Solves the puzzle by finding the transformation rule from the examples
    and applying it to the target input.
    """
    input_c = '000000,011120,111111'

    # Step 1: Parse the input string into a 2D grid of integers.
    rows_str = input_c.split(',')
    try:
        grid = [[int(char) for char in row] for row in rows_str]
    except ValueError:
        print("MISSING")
        return

    num_rows = len(grid)
    if num_rows == 0:
        print(input_c)
        return
    num_cols = len(grid[0])
    
    # Create a copy of the grid to store the result.
    # All transformation logic will read from the original `grid`.
    new_grid = [row[:] for row in grid]

    # Step 2: Find all locations of the number 2.
    twos_locations = []
    for r in range(num_rows):
        for c in range(num_cols):
            if grid[r][c] == 2:
                twos_locations.append((r, c))

    # If there are no 2s, the grid doesn't change.
    if not twos_locations:
        print(input_c)
        return

    # Step 3: For each '2', determine the transformation.
    for r, c in twos_locations:
        # Check if the '2' is at the left/right edge of the grid.
        if c == 0 or c == num_cols - 1:
            print("MISSING")
            return

        # Get Left and Right neighbors.
        L = grid[r][c - 1]
        R = grid[r][c + 1]

        # Determine jump direction based on the (L, R) pattern.
        dr, dc = None, None
        if (L, R) == (1, 1):
            dr, dc = -1, -1  # Top-left
        elif (L, R) == (0, 1):
            dr, dc = -1, 1   # Top-right
        elif (L, R) == (1, 0):
            dr, dc = 1, -1   # Bottom-left
        elif (L, R) == (0, 0):
            dr, dc = 1, 1    # Bottom-right
        else:
            # This handles cases where neighbors are not 0 or 1.
            print("MISSING")
            return
            
        # Calculate target coordinates.
        target_r, target_c = r + dr, c + dc

        # Validate the target: must be within bounds and its value must be 0.
        if not (0 <= target_r < num_rows and 0 <= target_c < num_cols):
            print("MISSING") # Target is out of bounds.
            return
        
        if grid[target_r][target_c] != 0:
            print("MISSING") # Target cell is not a '0'.
            return

        # Apply the transformation to the new grid.
        new_grid[r][c] = 0
        new_grid[target_r][target_c] = 2

    # If all transformations were valid, format and print the final grid.
    output_rows = ["".join(map(str, row)) for row in new_grid]
    final_output = ",".join(output_rows)
    print(final_output)

solve_puzzle()