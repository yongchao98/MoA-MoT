import numpy as np

def solve(grid):
    """
    Applies the transformation rule to the input grid.
    The new value of a cell is the sum of its 3x3 Moore neighborhood modulo 2.
    """
    rows, cols = grid.shape
    # Pad the grid with 0s to handle boundaries easily
    padded_grid = np.pad(grid, pad_width=1, mode='constant', constant_values=0)
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            # The 3x3 neighborhood in the padded grid corresponds to the
            # original cell's neighborhood.
            # The top-left corner of the subgrid is at (r, c) in the padded grid.
            neighborhood_sum = np.sum(padded_grid[r:r+3, c:c+3])
            output_grid[r, c] = neighborhood_sum % 2
    
    # Print the resulting grid as a flat string
    result_string = "".join(map(str, output_grid.flatten()))
    print(result_string)

# Test Input Grid 6
input_grid_6 = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

solve(input_grid_6)