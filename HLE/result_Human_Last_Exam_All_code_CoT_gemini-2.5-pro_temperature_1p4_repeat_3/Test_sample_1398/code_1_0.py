import numpy as np

def solve():
    """
    Solves the grid transformation puzzle by applying a specific rule.

    The rule is:
    A cell's value in the output grid is determined by the sum of its 8 neighbors in the input grid.
    - If the sum of the neighbors is 2 or 3, the output cell's value is 1.
    - Otherwise, the output cell's value is 0.
    The original value of the cell itself is ignored.
    """
    input_grid = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])

    rows, cols = input_grid.shape
    output_grid = np.zeros_like(input_grid)

    # Iterate through each cell of the grid
    for r in range(rows):
        for c in range(cols):
            # Define the boundaries for the neighborhood slice
            r_start = max(0, r - 1)
            r_end = min(rows, r + 2)
            c_start = max(0, c - 1)
            c_end = min(cols, c + 2)

            # Extract the 3x3 neighborhood
            neighborhood = input_grid[r_start:r_end, c_start:c_end]

            # Sum the neighbors (total sum of the neighborhood minus the cell itself)
            neighbor_sum = np.sum(neighborhood) - input_grid[r, c]

            # Apply the transformation rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    # Format the output grid as a single binary string
    output_string = "".join(map(str, output_grid.flatten()))
    
    print(output_string)

solve()