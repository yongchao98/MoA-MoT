def extract_subgrid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, input_grid[i:i+10])) for i in range(0, len(input_grid), 10)]
    
    # Define the size of the output grid
    output_size = 5
    
    # Find the center of the unique pattern
    # In this case, we will look for the distinct numbers and their positions
    distinct_positions = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != 3:  # Assuming 3 is the repeating number in this test case
                distinct_positions.append((i, j))
    
    # Calculate the center of the distinct positions
    if not distinct_positions:
        return []
    
    # Find the bounding box of the distinct positions
    min_i = min(pos[0] for pos in distinct_positions)
    max_i = max(pos[0] for pos in distinct_positions)
    min_j = min(pos[1] for pos in distinct_positions)
    max_j = max(pos[1] for pos in distinct_positions)
    
    # Calculate the center of the bounding box
    center_i = (min_i + max_i) // 2
    center_j = (min_j + max_j) // 2
    
    # Extract the subgrid centered around the calculated center
    start_i = max(0, center_i - output_size // 2)
    start_j = max(0, center_j - output_size // 2)
    
    output_grid = []
    for i in range(start_i, start_i + output_size):
        row = []
        for j in range(start_j, start_j + output_size):
            row.append(grid[i][j])
        output_grid.append(row)
    
    return output_grid

# Test input grid
input_grid = [
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 5, 3, 3, 5, 3, 3,
    3, 3, 3, 3, 3, 3, 4, 3, 3, 3,
    3, 3, 3, 3, 3, 9, 3, 3, 3, 3,
    3, 3, 0, 4, 5, 3, 3, 5, 3, 3,
    3, 3, 9, 0, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3
]

# Extract and print the output grid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))