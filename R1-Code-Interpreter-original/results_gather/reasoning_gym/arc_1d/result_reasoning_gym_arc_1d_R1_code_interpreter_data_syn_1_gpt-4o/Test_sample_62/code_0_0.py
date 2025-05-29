def transform_grid(input_grid):
    output_grid = []
    n = len(input_grid)
    for i in range(n):
        output_grid.append(input_grid[i])
        # Check if the current number is the end of a sequence of non-zero numbers
        if input_grid[i] != 0:
            # Check if the next number is zero or if we are at the end of the grid
            if i + 1 < n and input_grid[i + 1] == 0:
                output_grid.append(3)
            elif i + 1 == n:
                output_grid.append(3)
    return output_grid

# Test input
test_input = [4, 0, 0, 5, 5, 0, 0, 0, 9, 9, 9, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)