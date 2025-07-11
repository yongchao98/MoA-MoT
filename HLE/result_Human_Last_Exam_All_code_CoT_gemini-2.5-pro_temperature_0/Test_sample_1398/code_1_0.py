def solve_grid_puzzle():
    """
    Solves the grid transformation puzzle based on the neighbor sum rule.
    A cell in the output grid is 1 if the sum of its 8 neighbors in the
    input grid is 2 or 3. Otherwise, it's 0.
    """
    input_grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    height = len(input_grid)
    width = len(input_grid[0])
    output_grid = [[0 for _ in range(width)] for _ in range(height)]

    for r in range(height):
        for c in range(width):
            neighbor_sum = 0
            # Iterate through the 3x3 area centered on the cell
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    # Skip the cell itself
                    if dr == 0 and dc == 0:
                        continue

                    nr, nc = r + dr, c + dc

                    # Check if the neighbor is within the grid boundaries
                    if 0 <= nr < height and 0 <= nc < width:
                        neighbor_sum += input_grid[nr][nc]

            # Apply the transformation rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Print the final output grid row by row
    for row in output_grid:
        print(' '.join(map(str, row)))

solve_grid_puzzle()