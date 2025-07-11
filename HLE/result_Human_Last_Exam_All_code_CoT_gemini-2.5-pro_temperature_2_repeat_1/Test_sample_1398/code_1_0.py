import numpy as np

def solve_grid_puzzle():
    """
    Solves the grid puzzle by applying a transformation rule based on neighbor sums.
    The rule: a cell's output value is 1 if the sum of its 8 neighbors in the input grid is 2 or 3, otherwise it's 0.
    """
    
    # The 6th test input grid
    input_grid = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])

    rows, cols = input_grid.shape
    output_grid = np.zeros_like(input_grid)

    # Iterate over each cell to calculate its output value
    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of its 8 neighbors
            neighbor_sum = 0
            for i in range(max(0, r - 1), min(rows, r + 2)):
                for j in range(max(0, c - 1), min(cols, c + 2)):
                    # Skip the cell itself
                    if i == r and j == c:
                        continue
                    neighbor_sum += input_grid[i, j]

            # Apply the transformation rule
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    # The prompt requests to output each number in the final equation.
    # We will print the final resulting grid values, flattened into a single string.
    
    # Flatten the 2D grid into a 1D list of numbers
    final_numbers = output_grid.flatten()
    
    # Join the numbers into a single binary string
    binary_string_output = "".join(map(str, final_numbers))
    
    print(binary_string_output)

# Execute the function to get the solution
solve_grid_puzzle()