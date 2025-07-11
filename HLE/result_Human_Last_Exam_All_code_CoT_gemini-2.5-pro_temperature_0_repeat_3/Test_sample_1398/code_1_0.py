import collections

def solve():
    """
    This function solves the grid transformation puzzle.
    The rule, derived from the examples, is as follows:
    A cell in the output grid is 1 if the sum of its 8 neighbors in the input grid is 2 or 3.
    Otherwise, the cell is 0. The cell's own original value is irrelevant.
    This function applies this rule to the test grid and prints the result.
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

    for r in range(rows):
        for c in range(cols):
            neighbor_sum = 0
            # Iterate over the 8 neighbors
            for i in range(r - 1, r + 2):
                for j in range(c - 1, c + 2):
                    # Skip the cell itself
                    if i == r and j == c:
                        continue
                    # Check for boundary conditions
                    if 0 <= i < rows and 0 <= j < cols:
                        neighbor_sum += input_grid[i][j]
            
            # Apply the rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Print the output grid row by row
    for row in output_grid:
        print(" ".join(map(str, row)))

solve()