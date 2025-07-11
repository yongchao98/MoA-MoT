def solve_grid_transformation():
    """
    Solves the grid transformation problem by applying the discovered rule.
    The rule: A cell in the output grid becomes 1 if the sum of its 8 neighbors
    in the input grid is 2 or 3. Otherwise, it becomes 0. The original value
    of the cell is ignored.
    """
    
    # The 6th test input grid
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

    # Iterate through each cell of the grid to calculate its new state
    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of the 8 neighbors
            neighbor_sum = 0
            # Iterate over the 3x3 square centered at (r, c)
            for i in range(r - 1, r + 2):
                for j in range(c - 1, c + 2):
                    # Skip the cell itself
                    if i == r and j == c:
                        continue
                    # Check for boundary conditions before adding to the sum
                    if 0 <= i < rows and 0 <= j < cols:
                        neighbor_sum += input_grid[i][j]
            
            # Apply the transformation rule to determine the output cell's value
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Flatten the 2D output grid into a single string and print it.
    # This fulfills the request to "output each number in the final equation"
    # by printing all numbers of the resulting grid.
    output_string = ""
    for row in output_grid:
        for cell in row:
            output_string += str(cell)
    
    print(output_string)

# Execute the function to find the solution for the test grid
solve_grid_transformation()