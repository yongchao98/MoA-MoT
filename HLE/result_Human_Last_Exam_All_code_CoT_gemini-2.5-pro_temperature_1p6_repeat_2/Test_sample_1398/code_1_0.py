import numpy as np

def solve_grid_puzzle(grid):
    """
    Solves the grid puzzle based on the "B23/S" rule (a cell is 1 if it has 2 or 3 neighbors).
    This rule is a variation of Conway's Game of Life.
    """
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            # Count live neighbors in the 3x3 area
            live_neighbors = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    # Skip the cell itself
                    if dr == 0 and dc == 0:
                        continue
                    
                    nr, nc = r + dr, c + dc
                    
                    # Check if the neighbor is within the grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        live_neighbors += grid[nr, nc]
            
            # Apply the rule: a cell is 1 if it has 2 or 3 neighbors
            if live_neighbors == 2 or live_neighbors == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    return output_grid

# The 6th test input grid
input_grid_6 = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Calculate the output grid
output_grid_6 = solve_grid_puzzle(input_grid_6)

# Print the final grid as a 25-character binary string
flat_output = "".join(map(str, output_grid_6.flatten()))

print("The final grid is:")
for row in output_grid_6:
    print(" ".join(map(str, row)))

# The problem asks for the output in a specific format in the final response,
# but for the code block, it's helpful to also see the grid structure.
# Here we will just print the final flattened string as requested by the user prompt
# to make sure the format is correct for the final answer submission.
# print("\nFinal binary string:")
# print(flat_output)
