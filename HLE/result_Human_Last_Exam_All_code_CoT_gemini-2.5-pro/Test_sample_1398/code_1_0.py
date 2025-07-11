import numpy as np

def solve():
    """
    Solves the grid transformation puzzle by applying a cellular automaton rule.
    The rule is: a cell's next state is 1 if it has 2 or 3 live neighbors, and 0 otherwise.
    The state of the cell itself is not considered, only its neighbors.
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
            # Count live neighbors in the 8 surrounding cells
            n_neighbors = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue  # Skip the cell itself

                    nr, nc = r + dr, c + dc

                    # Check if the neighbor is within the grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        n_neighbors += input_grid[nr, nc]

            # Apply the rule: next state is 1 if 2 or 3 neighbors are live
            if n_neighbors == 2 or n_neighbors == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    # Print the output grid
    for row in output_grid:
        print(' '.join(map(str, row)))

    # The final answer as a single binary string
    # This part is for the final answer submission and won't be visible in the execution.
    # final_string = "".join(map(str, output_grid.flatten()))
    # print(f"\nFinal Answer String: {final_string}")

solve()