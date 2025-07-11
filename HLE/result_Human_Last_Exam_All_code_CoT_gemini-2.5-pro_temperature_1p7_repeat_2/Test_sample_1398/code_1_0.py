import numpy as np

def solve_grid(grid):
    """
    Applies the B235/S23 cellular automaton rule to a grid.
    - A dead cell (0) becomes live if it has 2, 3, or 5 live neighbors.
    - A live cell (1) survives if it has 2 or 3 live neighbors.
    """
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)
    
    # Define birth and survival conditions
    birth_conditions = {2, 3, 5}
    survival_conditions = {2, 3}
    
    # Iterate through each cell in the grid
    for r in range(rows):
        for c in range(cols):
            # Sum of the 8 neighbors (Moore neighborhood)
            live_neighbors = 0
            for i in range(-1, 2):
                for j in range(-1, 2):
                    # Skip the cell itself
                    if i == 0 and j == 0:
                        continue
                    
                    nr, nc = r + i, c + j
                    
                    # Check grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        live_neighbors += grid[nr, nc]

            # Apply the B235/S23 rule
            current_cell_value = grid[r, c]
            if current_cell_value == 0: # Birth rule
                if live_neighbors in birth_conditions:
                    output_grid[r, c] = 1
            elif current_cell_value == 1: # Survival rule
                if live_neighbors in survival_conditions:
                    output_grid[r, c] = 1
            # Otherwise, the cell is 0 by default initialization

    return output_grid

def print_grid(grid):
    """Prints the grid in the required format."""
    for row in grid:
        # The prompt mentions "output each number in the final equation"
        # which is interpreted as printing the grid values clearly.
        print(' '.join(map(str, row)))

# Test Input Grid 6
input_grid_6 = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Calculate and print the output grid
output_grid_6 = solve_grid(input_grid_6)
print_grid(output_grid_6)

# The following line is for generating the final answer string, not part of the script's visible output.
# final_answer_string = "".join(map(str, output_grid_6.flatten()))
# print(f"<<<{final_answer_string}>>>")