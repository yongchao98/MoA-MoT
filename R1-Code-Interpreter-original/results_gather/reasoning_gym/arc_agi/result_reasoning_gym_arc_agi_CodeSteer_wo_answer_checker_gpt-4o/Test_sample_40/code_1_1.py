def transform_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    rows = len(input_grid)
    cols = len(input_grid[0])

    # Directions for adjacency (up, down, left, right)
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    # Function to check if a position is within the grid
    def in_bounds(x, y):
        return 0 <= x < rows and 0 <= y < cols

    # Iterate over each cell in the grid
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 0:
                # Check for specific patterns around the cell
                adjacent_fours = 0
                for dx, dy in directions:
                    ni, nj = i + dx, j + dy
                    if in_bounds(ni, nj) and input_grid[ni][nj] == 4:
                        adjacent_fours += 1

                # Apply pattern-based rules
                if adjacent_fours >= 2:
                    # If there are two or more adjacent '4's, place a '3'
                    output_grid[i][j] = 3
                elif adjacent_fours == 1:
                    # If there is exactly one adjacent '4', place a '1'
                    output_grid[i][j] = 1

    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0],
    [0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 6, 4, 4, 0, 0],
    [0, 0, 4, 0, 4, 0, 0, 0, 0, 1, 4, 0, 4, 0, 0],
    [0, 0, 6, 4, 4, 0, 0, 0, 0, 3, 4, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 4, 4, 6, 0, 0, 0, 0, 0, 4, 4, 0],
    [0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 4, 0],
    [0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 4, 4, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))