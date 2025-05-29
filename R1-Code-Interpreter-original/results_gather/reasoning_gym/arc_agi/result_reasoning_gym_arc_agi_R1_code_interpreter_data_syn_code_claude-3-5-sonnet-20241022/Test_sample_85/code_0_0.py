def find_positions(grid, target):
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] == target:
                return (i, j)
    return None

def create_path(grid, start, end):
    result = [row[:] for row in grid]
    
    # Get coordinates
    start_row, start_col = start
    end_row, end_col = end
    
    # Create horizontal path
    if start_row == end_row:
        for col in range(min(start_col, end_col) + 1, max(start_col, end_col)):
            result[start_row][col] = 8
    else:
        # Create vertical path from start to end's column
        current_row = start_row
        direction = 1 if end_row > start_row else -1
        while current_row != end_row:
            result[current_row][start_col] = 8
            current_row += direction
        
        # Create horizontal path
        for col in range(min(start_col, end_col), max(start_col, end_col) + 1):
            if result[end_row][col] not in [1, 9]:
                result[end_row][col] = 8
    
    return result

# Test input
test_input = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Find positions of 1 and 9
pos_1 = find_positions(test_input, 1)
pos_9 = find_positions(test_input, 9)

# Create output grid
output = create_path(test_input, pos_1, pos_9)

# Print the result in the required format
print("<<<")
for row in output:
    print(" ".join(map(str, row)))
print(">>>")