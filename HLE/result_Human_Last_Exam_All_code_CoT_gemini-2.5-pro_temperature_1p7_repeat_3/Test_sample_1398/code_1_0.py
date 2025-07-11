import numpy as np

def solve_grid(grid):
    """
    Applies the HighLife (B23/S23) cellular automaton rule to a grid.
    - A dead cell (0) becomes alive if it has 2 or 3 neighbors (3x3 sum is 2 or 3).
    - An alive cell (1) survives if it has 2 or 3 neighbors (3x3 sum is 3 or 4).
    """
    # Pad the grid with zeros to handle boundaries
    padded_grid = np.pad(grid, pad_width=1, mode='constant', constant_values=0)
    output_grid = np.zeros_like(grid)
    
    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            # Sum the 3x3 neighborhood in the padded grid
            # The coordinates in the padded grid are (r+1, c+1)
            neighborhood_sum = np.sum(padded_grid[r:r+3, c:c+3])
            
            current_val = grid[r, c]
            
            # Apply the B23/S23 rule
            if current_val == 0:  # If cell is dead
                if neighborhood_sum == 2 or neighborhood_sum == 3:
                    output_grid[r, c] = 1
            else:  # If cell is alive
                if neighborhood_sum == 3 or neighborhood_sum == 4:
                    output_grid[r, c] = 1
                    
    return output_grid

def main():
    # Test Input Grid 6
    input_grid = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])
    
    output_grid = solve_grid(input_grid)
    
    # Print the output grid
    for row in output_grid:
        print(' '.join(map(str, row)))

main()
