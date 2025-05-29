def find_patterns(grid):
    patterns = []
    rows, cols = len(grid), len(grid[0])
    
    # Find horizontal patterns
    for i in range(rows):
        j = 0
        while j < cols:
            if grid[i][j] != 4:
                # Find length of pattern
                start = j
                current_value = grid[i][j]
                pattern = []
                while j < cols and grid[i][j] != 4 and grid[i][j] == current_value:
                    pattern.append(grid[i][j])
                    j += 1
                if len(pattern) > 1:  # Only consider patterns of length > 1
                    patterns.append(('h', i, start, j-1, pattern))
            j += 1
    
    return patterns

def transform_pattern(pattern):
    type, row, start, end, values = pattern
    length = end - start + 1
    
    # Only transform patterns of 8s
    if all(v == 8 for v in values):
        if length == 2:
            return [0, 0]
        elif length == 3:
            return [9, 9, 9]
        elif length == 4:
            return [9, 9, 9, 9]
    
    return values  # If no transformation rule matches, return original values

def solve_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy input grid
    
    patterns = find_patterns(input_grid)
    
    for pattern in patterns:
        type, row, start, end, values = pattern
        transformed = transform_pattern(pattern)
        for i, val in enumerate(transformed):
            output_grid[row][start + i] = val
    
    # Print the result in the required format
    result = ""
    for row in output_grid:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
input_grid = [
    [4,4,4,4,4,4,4,4,4,4],
    [4,4,8,8,4,4,4,8,8,4],
    [4,8,8,8,8,4,8,8,8,8],
    [4,4,8,4,4,4,4,8,4,4],
    [4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,2,2,4],
    [4,4,8,8,4,4,1,1,2,1],
    [4,8,8,8,8,4,4,1,4,4],
    [4,4,8,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4]
]

solve_grid(input_grid)