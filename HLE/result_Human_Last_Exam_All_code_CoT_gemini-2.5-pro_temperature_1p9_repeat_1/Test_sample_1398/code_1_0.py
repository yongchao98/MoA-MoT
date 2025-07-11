def solve_grid_puzzle(grid):
    """
    Applies the transformation rule to a grid.
    The rule is: A cell's new state is 1 if it has 2 or 3 live neighbors,
    and 0 otherwise. The original state of the cell does not matter.
    """
    if not grid or not grid[0]:
        return ""

    rows = len(grid)
    cols = len(grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    for r in range(rows):
        for c in range(cols):
            live_neighbors = 0
            # Iterate through the 8 neighbors
            for i in range(r - 1, r + 2):
                for j in range(c - 1, c + 2):
                    # Skip the cell itself
                    if i == r and j == c:
                        continue
                    # Check for boundary conditions
                    if 0 <= i < rows and 0 <= j < cols:
                        if grid[i][j] == 1:
                            live_neighbors += 1
            
            # The final equation is: new_value = 1 if live_neighbors is 2 or 3, else 0.
            # The numbers in the equation are 2 and 3.
            if live_neighbors == 2 or live_neighbors == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Flatten the grid into a single binary string
    result_string = "".join(str(cell) for row in output_grid for cell in row)
    print(result_string)

# Test Input Grid 6
test_input_grid = [
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
]

solve_grid_puzzle(test_input_grid)