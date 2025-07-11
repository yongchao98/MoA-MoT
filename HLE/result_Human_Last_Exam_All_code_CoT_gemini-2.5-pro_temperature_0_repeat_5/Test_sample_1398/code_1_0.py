import numpy as np

def solve_grid():
    """
    Solves the grid transformation puzzle by applying a discovered rule.
    The rule is a variant of Conway's Game of Life: B2356/S23.
    - A dead cell (0) becomes live (1) if it has 2, 3, 5, or 6 live neighbors.
    - A live cell (1) survives (stays 1) if it has 2 or 3 live neighbors.
    - Otherwise, a cell is dead (0) in the next state.
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

    for r in range(rows):
        for c in range(cols):
            # Define the 3x3 neighborhood boundaries, handling edges
            r_start, r_end = max(0, r - 1), min(rows, r + 2)
            c_start, c_end = max(0, c - 1), min(cols, c + 2)
            
            # Sum the neighbors
            neighbor_sum = np.sum(input_grid[r_start:r_end, c_start:c_end])
            
            # The sum includes the cell itself, so subtract it
            center_val = input_grid[r, c]
            neighbor_sum -= center_val

            # Apply the B2356/S23 rule
            if center_val == 0:  # Dead cell (Birth rule)
                if neighbor_sum in [2, 3, 5, 6]:
                    output_grid[r, c] = 1
            else:  # Live cell (Survival rule)
                if neighbor_sum in [2, 3]:
                    output_grid[r, c] = 1
    
    # Print the final grid as requested
    for row in output_grid:
        print(' '.join(map(str, row)))

solve_grid()