def transform_grid(input_grid):
    # Determine the most common number in the grid
    from collections import Counter
    flat_grid = [num for row in input_grid for num in row]
    most_common_num = Counter(flat_grid).most_common(1)[0][0]
    
    # Create a copy of the input grid for the output
    output_grid = [row[:] for row in input_grid]
    
    # Define the transformation rules based on the examples
    transformation_rules = {
        7: 7,
        9: 4,
        6: 5,
        8: 5,
        5: 9
    }
    
    # Apply the transformation
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if input_grid[i][j] in transformation_rules:
                # Get the transformation number
                transform_num = transformation_rules[input_grid[i][j]]
                # Surround the number with the transform_num
                for di in [-1, 0, 1]:
                    for dj in [-1, 0, 1]:
                        ni, nj = i + di, j + dj
                        if 0 <= ni < len(input_grid) and 0 <= nj < len(input_grid[0]):
                            if input_grid[ni][nj] == most_common_num:
                                output_grid[ni][nj] = transform_num
    
    return output_grid

# Test input grid
input_grid = [
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 7, 2, 2, 2, 2, 7, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [7, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 7, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))