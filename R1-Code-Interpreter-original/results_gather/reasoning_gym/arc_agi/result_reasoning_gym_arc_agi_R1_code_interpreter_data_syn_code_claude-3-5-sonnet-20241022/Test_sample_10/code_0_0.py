def find_patterns(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Function to check if a position contains a 3x3 pattern with a 9
    def has_3x3_pattern(r, c):
        if r+2 >= rows or c+2 >= cols:
            return False
        count_1s = 0
        has_9 = False
        for i in range(r, r+3):
            for j in range(c, c+3):
                if input_grid[i][j] == 1:
                    count_1s += 1
                elif input_grid[i][j] == 9:
                    has_9 = True
        return count_1s >= 7 and has_9

    # Find all 3x3 patterns in input
    patterns = []
    for i in range(rows):
        for j in range(cols):
            if has_3x3_pattern(i, j):
                patterns.append((i, j))
                
    # Check which patterns remain in output
    remaining_patterns = []
    for i, j in patterns:
        if has_3x3_pattern(i, j):
            if all(output_grid[r][c] == input_grid[r][c] 
                  for r in range(i, i+3) 
                  for c in range(j, j+3)):
                remaining_patterns.append((i, j))
    
    print(f"Total patterns found in input: {len(patterns)}")
    print(f"Patterns remaining in output: {len(remaining_patterns)}")
    return patterns, remaining_patterns

# Example input and output from Example 1
input1 = [
    [8,8,1,1,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,9,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,1,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,1,1,8,8,1,1,1,1,8,8,8],
    [8,8,8,8,8,8,8,1,1,9,1,8,8,8],
    [1,1,1,8,8,8,8,1,9,1,1,8,8,8],
    [9,1,1,8,8,8,8,1,1,1,1,8,8,8],
    [1,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,1,1,1,8,8,8,8,8,8],
    [8,8,8,8,8,1,9,1,8,8,1,1,1,1],
    [1,1,1,8,8,9,1,1,8,8,9,1,9,1],
    [1,1,1,8,8,1,1,9,8,8,1,1,9,1],
    [8,8,8,8,8,8,8,8,8,8,1,9,1,1]
]

output1 = [
    [8,8,1,1,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,9,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,1,1,8,8,8,8,8,8,8,8,8],
    [8,8,1,1,1,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8],
    [1,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [9,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [1,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8],
    [1,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [1,1,1,8,8,8,8,8,8,8,8,8,8,8],
    [8,8,8,8,8,8,8,8,8,8,8,8,8,8]
]

find_patterns(input1, output1)