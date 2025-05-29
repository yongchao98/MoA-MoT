def transform_grid(input_grid):
    # Initialize the output grid with zeros
    output_grid = [[0] * 10 for _ in range(11)]
    
    # Place '8' at the corners
    output_grid[0][0] = 8
    output_grid[0][9] = 8
    output_grid[10][0] = 8
    output_grid[10][9] = 8
    
    # Identify unique numbers in the input grid
    unique_numbers = set()
    for row in input_grid:
        unique_numbers.update(row)
    
    # Remove '0' and '8' from unique numbers
    unique_numbers.discard(0)
    unique_numbers.discard(8)
    
    # Sort the remaining numbers to find the largest and second largest
    sorted_numbers = sorted(unique_numbers, reverse=True)
    
    if len(sorted_numbers) >= 1:
        max_num = sorted_numbers[0]
    else:
        max_num = 0
    
    if len(sorted_numbers) >= 2:
        second_max_num = sorted_numbers[1]
    else:
        second_max_num = 0
    
    # Place the largest number in the middle block
    for i in range(3, 7):
        for j in range(3, 7):
            output_grid[i][j] = max_num
    
    # Place the second largest number in the lower part
    for i in range(7, 10):
        for j in range(3, 7):
            output_grid[i][j] = second_max_num
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))