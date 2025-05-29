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
    
    # First create horizontal path from 1 to 9's column
    for col in range(min(start_col + 1, end_col), max(start_col, end_col) + 1):
        if result[start_row][col] != 9:  # Don't overwrite 9
            result[start_row][col] = 8
            
    # Then create vertical path to 9
    if end_row < start_row:  # If 9 is above
        for row in range(end_row, start_row):
            if result[row][end_col] != 9:  # Don't overwrite 9
                result[row][end_col] = 8
    else:  # If 9 is below
        for row in range(start_row + 1, end_row + 1):
            if result[row][end_col] != 9:  # Don't overwrite 9
                result[row][end_col] = 8
    
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