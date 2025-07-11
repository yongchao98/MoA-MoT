import collections

def solve_grid():
    """
    Solves the grid transformation puzzle by discovering the rule from examples
    and applying it to the test grid.
    """

    test_grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    def get_neighbor_sum(grid, r, c):
        """Calculates the sum of the 8 neighbors for a cell at (r, c)."""
        rows, cols = len(grid), len(grid[0])
        s = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    s += grid[nr][nc]
        return s

    # The rule is derived from analyzing the provided examples.
    # Output is 1 if (input_value, neighbor_sum) matches specific pairs.
    # Specifically: O=1 if (I=0, N=2), (I=0, N=3), or (I=1, N=2).
    # Otherwise, O=0.
    
    rows, cols = len(test_grid), len(test_grid[0])
    output_grid = [[0] * cols for _ in range(rows)]

    for r in range(rows):
        for c in range(cols):
            n_sum = get_neighbor_sum(test_grid, r, c)
            input_val = test_grid[r][c]

            # Apply the discovered transformation rule
            if (input_val == 0 and n_sum in [2, 3]) or \
               (input_val == 1 and n_sum == 2):
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Format the output as a single 25-character binary string
    final_string = ""
    for row in output_grid:
        final_string += "".join(map(str, row))

    print(final_string)

solve_grid()