import numpy as np

def solve():
    """
    Solves the grid transformation puzzle.
    The transformation rule is based on a cell's own value (C) and the sum of its
    cardinal (S4) and diagonal (Sd) neighbors.
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
    
    # Pad the grid with 0s to handle boundaries easily
    padded_grid = np.pad(input_grid, pad_width=1, mode='constant', constant_values=0)

    for r in range(rows):
        for c in range(cols):
            # Coordinates in the padded grid
            pr, pc = r + 1, c + 1
            
            # C: Value of the current cell
            C = padded_grid[pr, pc]

            # S4: Sum of cardinal neighbors
            S4 = padded_grid[pr-1, pc] + padded_grid[pr+1, pc] + \
                 padded_grid[pr, pc-1] + padded_grid[pr, pc+1]
            
            # Sd: Sum of diagonal neighbors
            Sd = padded_grid[pr-1, pc-1] + padded_grid[pr-1, pc+1] + \
                 padded_grid[pr+1, pc-1] + padded_grid[pr+1, pc+1]

            output_val = 0
            if C == 0:
                if S4 == 2 and Sd == 0:
                    output_val = 1
                else:
                    output_val = 0
            elif C == 1:
                if (S4, Sd) in [(0, 2), (1, 1), (1, 3), (2, 2)]:
                    output_val = 1
                else:
                    output_val = 0
            
            output_grid[r, c] = output_val
    
    # Flatten the grid and print as a binary string
    binary_string = "".join(map(str, output_grid.flatten()))
    print(f"The final binary string is: {binary_string}")

solve()