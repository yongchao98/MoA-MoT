def solve_puzzle():
    """
    Solves the puzzle by finding and applying the transformation rule.
    The rule involves moving the digit '2' based on its neighboring '1's.
    """
    input_c = '000000,011120,111111'

    # Parse the input string into a 2D list of integers.
    grid_str_list = input_c.split(',')
    grid = [[int(char) for char in row] for row in grid_str_list]
    height = len(grid)
    width = len(grid[0])

    # Find the coordinates (r, c) of the '2'.
    r2, c2 = -1, -1
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    # Define the 8 neighbor offsets.
    neighbor_offsets = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),          (0, 1),
        (1, -1), (1, 0), (1, 1)
    ]

    # Count the number of '1's in the neighborhood of '2'.
    num_ones = 0
    for dr, dc in neighbor_offsets:
        nr, nc = r2 + dr, c2 + dc
        if 0 <= nr < height and 0 <= nc < width:
            if grid[nr][nc] == 1:
                num_ones += 1

    # Define the direction mapping based on the discovered rule.
    # The order is a top-to-bottom, right-to-left scan of neighbors.
    directions = [
        (-1, 1), (-1, 0), (-1, -1),
        (0, 1),           (0, -1),
        (1, 1), (1, 0), (1, -1)
    ]
    
    # Calculate the direction index based on the number of '1's.
    direction_index = (num_ones - 4) % 8
    
    # Get the target direction and calculate the new coordinates.
    dr_target, dc_target = directions[direction_index]
    tr, tc = r2 + dr_target, c2 + dc_target
    
    # Apply the transformation to the grid.
    grid[r2][c2] = 0
    grid[tr][tc] = 2
    
    # Convert the modified grid back to the comma-separated string format.
    result_rows = ["".join(map(str, row)) for row in grid]
    final_output = ",".join(result_rows)

    print(final_output)

solve_puzzle()