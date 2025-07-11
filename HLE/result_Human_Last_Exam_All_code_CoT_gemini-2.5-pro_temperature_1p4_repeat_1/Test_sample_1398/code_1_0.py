import numpy as np

def solve_grid_puzzle(grid):
    """
    Solves the grid puzzle based on the neighbor counting rule.
    A cell's new state is 1 if it has 2 or 3 live neighbors, and 0 otherwise.
    """
    # Get the dimensions of the grid
    rows, cols = grid.shape
    
    # Create a new grid to store the output, initialized with zeros
    output_grid = np.zeros_like(grid)
    
    # Iterate over each cell in the input grid
    for r in range(rows):
        for c in range(cols):
            # Count the number of live neighbors (value 1)
            live_neighbors = 0
            # Iterate over the 3x3 neighborhood around the cell
            for i in range(max(0, r - 1), min(rows, r + 2)):
                for j in range(max(0, c - 1), min(cols, c + 2)):
                    # Skip the cell itself
                    if i == r and j == c:
                        continue
                    # If the neighbor is within bounds and is 'live', increment count
                    if grid[i, j] == 1:
                        live_neighbors += 1
            
            # Apply the rule: new cell is 1 if live_neighbors is 2 or 3
            if live_neighbors == 2 or live_neighbors == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0
                
    return output_grid

def print_grid(grid):
    """
    Prints the grid in the required format.
    """
    for row in grid:
        print(' '.join(map(str, row)))

# The 6th input grid
test_input_grid = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Solve the puzzle for the given grid
final_output_grid = solve_grid_puzzle(test_input_grid)

# Print the final grid
print_grid(final_output_grid)