def transform_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Initialize the output grid with 2s
    output_size = max(input_rows, input_cols) * 3
    output_grid = [[2 for _ in range(output_size)] for _ in range(output_size)]
    
    # Identify unique numbers in the input grid
    unique_numbers = set()
    for row in input_grid:
        unique_numbers.update(row)
    
    # Sort unique numbers in descending order to prioritize larger numbers
    sorted_numbers = sorted(unique_numbers, reverse=True)
    
    # Map each unique number to a block in the output grid
    block_size = 3
    for num in sorted_numbers:
        # Find the first occurrence of the number in the input grid
        for i in range(input_rows):
            for j in range(input_cols):
                if input_grid[i][j] == num:
                    # Place a block in the output grid
                    for x in range(block_size):
                        for y in range(block_size):
                            output_grid[i * block_size + x][j * block_size + y] = num
                    break
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))