import numpy as np

def solve():
    """
    Solves the grid transformation puzzle.
    The transformation rule is determined to be:
    The new value of a cell is the sum of its original value and the values
    of its four diagonal neighbors, all taken modulo 2.
    O(r, c) = (I(r, c) + sum(diagonal_neighbors)) % 2
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

    # Pad the grid with 0s to handle boundary conditions easily
    padded_grid = np.pad(input_grid, pad_width=1, mode='constant', constant_values=0)

    for r in range(rows):
        for c in range(cols):
            # Padded grid coordinates are (r+1, c+1)
            pr, pc = r + 1, c + 1

            # Sum of diagonal neighbors
            s_d = (padded_grid[pr - 1, pc - 1] +
                   padded_grid[pr - 1, pc + 1] +
                   padded_grid[pr + 1, pc - 1] +
                   padded_grid[pr + 1, pc + 1])

            # Current cell's value
            i_rc = padded_grid[pr, pc]

            # Apply the rule
            output_grid[r, c] = (i_rc + s_d) % 2
            
    # Print the output grid row by row
    # The requirement "output each number in the final equation" is interpreted
    # as printing the final grid which constitutes the answer.
    final_string = ""
    for r in range(rows):
        row_str = " ".join(map(str, output_grid[r]))
        print(row_str)
        final_string += "".join(map(str, output_grid[r]))

    print(f"\nFinal binary string:\n{final_string}")


solve()