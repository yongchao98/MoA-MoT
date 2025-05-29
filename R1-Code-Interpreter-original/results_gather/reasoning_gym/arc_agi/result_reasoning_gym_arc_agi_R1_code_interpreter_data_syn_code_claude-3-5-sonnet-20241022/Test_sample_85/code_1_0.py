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
    
    # First create vertical path from start to end's row
    current_row = min(start_row, end_row)
    end_vert = max(start_row, end_row)
    
    # If 9 is above 1, go up from 1
    if end_row < start_row:
        for row in range(start_row - 1, end_row - 1, -1):
            result[row][start_col] = 8
    # If 9 is below 1 or at same level, go up from 1's position
    else:
        for row in range(start_row + 1, end_row + 1):
            result[row][start_col] = 8
    
    # Then create horizontal path at end's row
    min_col = min(start_col, end_col)
    max_col = max(start_col, end_col)
    for col in range(min_col, max_col + 1):
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