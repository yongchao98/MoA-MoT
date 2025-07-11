import numpy as np

def solve():
    """
    This function solves the grid transformation puzzle.
    The transformation rule is that a cell's new state is 1 if the sum of its 8 neighbors 
    in the input grid is 2 or 3, and 0 otherwise.
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
            neighbor_sum = 0
            # Iterate through the 8 neighbors
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    # Skip the cell itself
                    if dr == 0 and dc == 0:
                        continue

                    nr, nc = r + dr, c + dc

                    # Check if the neighbor is within grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        neighbor_sum += input_grid[nr, nc]

            # Apply the transformation rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    # Print the resulting output grid
    for row in output_grid:
        print(' '.join(map(str, row)))

    # Generate and prepare the final answer string
    binary_string = "".join(map(str, output_grid.flatten()))
    # The prompt asks for a specific format for the final answer.
    # We print it here directly as requested.
    print(f"\n<<<binary_string:{binary_string}>>>")


solve()