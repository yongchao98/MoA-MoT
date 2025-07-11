import numpy as np

def solve_grid(grid):
    """
    Solves the grid transformation puzzle.
    The rule is: O(r, c) = I(r, c) XOR (S(r, c) % 2)
    where I is the input value and S is the sum of Von Neumann neighbors.
    """
    input_grid = np.array(grid)
    output_grid = np.zeros_like(input_grid)
    rows, cols = input_grid.shape

    for r in range(rows):
        for c in range(cols):
            neighbor_sum = 0
            # North neighbor
            if r > 0:
                neighbor_sum += input_grid[r - 1, c]
            # South neighbor
            if r < rows - 1:
                neighbor_sum += input_grid[r + 1, c]
            # West neighbor
            if c > 0:
                neighbor_sum += input_grid[r, c - 1]
            # East neighbor
            if c < cols - 1:
                neighbor_sum += input_grid[r, c + 1]

            # Apply the transformation rule
            parity = neighbor_sum % 2
            new_value = input_grid[r, c] ^ parity
            output_grid[r, c] = new_value

    # Print the output grid as a flat binary string
    flat_string = "".join(map(str, output_grid.flatten()))
    print(flat_string)

# The 6th test input grid
test_input_grid = [
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
]

solve_grid(test_input_grid)