def find_special_rows(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    # Count frequency of each row pattern
    row_patterns = {}
    for i in range(rows):
        row_tuple = tuple(grid[i])
        row_patterns[row_tuple] = row_patterns.get(row_tuple, 0) + 1
    
    # Find the majority pattern and special rows
    majority_pattern = max(row_patterns.items(), key=lambda x: x[1])[0]
    special_rows = []
    for i in range(rows):
        if tuple(grid[i]) != majority_pattern:
            special_rows.append(i)
    
    return special_rows

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    special_rows = find_special_rows(input_grid)
    
    # If we have consecutive special rows
    if len(special_rows) >= 2 and special_rows[1] == special_rows[0] + 1:
        # For each special row
        for row in special_rows:
            # For each column
            for col in range(cols):
                special_value = input_grid[row][col]  # Value in special row
                vertical_value = input_grid[0][col]  # Value from top row
                # If special value is larger than vertical value, inherit vertical value
                if special_value > vertical_value:
                    output_grid[row][col] = vertical_value
                # If special value is smaller or equal, keep original value
                else:
                    output_grid[row][col] = special_value
    
    # Print the result
    for row in output_grid:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [2,2,7,7,2,2,2,2,2,2,2],
    [2,2,7,7,2,2,2,2,2,2,2],
    [3,3,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3,3],
    [2,2,7,7,2,2,2,2,2,2,2],
    [2,2,7,7,2,2,2,2,2,2,2]
]

transform_grid(test_input)