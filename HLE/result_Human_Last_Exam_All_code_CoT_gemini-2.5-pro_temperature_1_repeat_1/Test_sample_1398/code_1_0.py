import numpy as np

def solve_grid_puzzle(grid):
    """
    Applies the transformation rule to a grid.
    The rule: a cell's output value is 1 if the sum of its 8 neighbors in the input is 2 or 3, else 0.
    """
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            # Define the boundaries for the neighborhood
            r_min = max(0, r - 1)
            r_max = min(rows, r + 2)
            c_min = max(0, c - 1)
            c_max = min(cols, c + 2)

            # Extract the 3x3 neighborhood
            neighborhood = grid[r_min:r_max, c_min:c_max]
            
            # Sum all values in the neighborhood
            neighbor_sum = np.sum(neighborhood)
            
            # Subtract the cell's own value from the sum
            neighbor_sum -= grid[r, c]

            # Apply the rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0
                
    return output_grid

# Test Input Grid 6
input_grid = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Calculate the output grid
output_grid = solve_grid_puzzle(input_grid)

# Print the output grid in a readable format
print("The calculated output grid is:")
for row in output_grid:
    # This loop prints each number in the final grid, as requested by the prompt ("output each number in the final equation").
    print(' '.join(map(str, row)))

# Generate the final binary string
binary_string = "".join(map(str, output_grid.flatten()))
print("\nThe final answer as a 25-value binary string:")
print(binary_string)