import numpy as np

def solve(grid):
    """
    Applies the transformation rule to solve the grid puzzle.
    """
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            # Get neighbor coordinates, handling boundaries by considering them 0
            # A helper function to safely get grid values
            def get_val(r_idx, c_idx):
                if 0 <= r_idx < rows and 0 <= c_idx < cols:
                    return grid[r_idx, c_idx]
                return 0

            # Get values of the 4 diagonal neighbors
            diag_neighbors = [
                get_val(r - 1, c - 1),
                get_val(r - 1, c + 1),
                get_val(r + 1, c - 1),
                get_val(r + 1, c + 1)
            ]
            sum_diag = sum(diag_neighbors)

            # Get values of the 2 horizontal neighbors
            h_left = get_val(r, c - 1)
            h_right = get_val(r, c + 1)
            
            current_cell_val = grid[r, c]

            # Apply the special horizontal rule first
            if current_cell_val == 1 and h_left == 1 and h_right == 1:
                output_grid[r, c] = 1
            # Otherwise, apply the default diagonal rule
            else:
                output_grid[r, c] = sum_diag % 2
    
    return output_grid

def main():
    # Test Input Grid 6
    input_grid_6 = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])

    # Solve for the output grid
    output_grid_6 = solve(input_grid_6)
    
    # Print the result as a flattened binary string
    result_string = "".join(map(str, output_grid_6.flatten()))
    print(result_string)

main()