def solve_grid_puzzle():
    """
    Solves the grid puzzle by applying the deduced transformation rule.
    The rule depends on the cell's value and the parity of its von Neumann and Moore neighbor sums.
    """
    input_grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    def get_val(grid, r, c):
        """Safely gets a value from the grid, returning 0 for out-of-bounds indices."""
        if 0 <= r < rows and 0 <= c < cols:
            return grid[r][c]
        return 0

    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of the four direct (von Neumann) neighbors
            n_sum = (get_val(input_grid, r - 1, c) +
                     get_val(input_grid, r + 1, c) +
                     get_val(input_grid, r, c - 1) +
                     get_val(input_grid, r, c + 1))

            # Calculate the sum of the eight surrounding (Moore) neighbors
            m_sum = (get_val(input_grid, r - 1, c - 1) + get_val(input_grid, r - 1, c) + get_val(input_grid, r - 1, c + 1) +
                     get_val(input_grid, r, c - 1) + get_val(input_grid, r, c + 1) +
                     get_val(input_grid, r + 1, c - 1) + get_val(input_grid, r + 1, c) + get_val(input_grid, r + 1, c + 1))

            current_val = get_val(input_grid, r, c)

            # Apply the deduced transformation rule
            if current_val == 0:
                output_grid[r][c] = m_sum % 2
            else:  # current_val == 1
                if n_sum % 2 == 0:
                    output_grid[r][c] = 0
                else:
                    output_grid[r][c] = 1 - (m_sum % 2)

    # Print the final result as a 25-character binary string
    binary_string = "".join(str(cell) for row in output_grid for cell in row)
    print(binary_string)

solve_grid_puzzle()