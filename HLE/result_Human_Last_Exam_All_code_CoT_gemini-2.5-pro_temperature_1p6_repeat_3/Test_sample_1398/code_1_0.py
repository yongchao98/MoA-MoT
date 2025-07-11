import numpy as np

def solve_grid(grid):
    """
    Applies the transformation rule to solve the grid puzzle.
    The rule: a cell's new value is 1 minus the sum of its 3x3 neighborhood modulo 2.
    """
    rows, cols = grid.shape
    # Pad the grid with a border of 0s to handle edge cases easily
    padded_grid = np.pad(grid, pad_width=1, mode='constant', constant_values=0)
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            # We use r+1 and c+1 because of the padding
            pr, pc = r + 1, c + 1
            
            # Sum the values in the 3x3 neighborhood in the padded grid
            neighborhood_sum = np.sum(padded_grid[pr-1:pr+2, pc-1:pc+2])
            
            # Apply the rule: 1 - (sum mod 2)
            output_grid[r, c] = 1 - (neighborhood_sum % 2)
            
    return output_grid

# The sixth input grid
test_input_grid = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Calculate the output grid
result_grid = solve_grid(test_input_grid)

# Print the final output grid as a flat binary string
flat_result = "".join(map(str, result_grid.flatten()))
print(flat_result)