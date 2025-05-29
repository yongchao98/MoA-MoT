def find_9_pattern(grid):
    # Find boundaries of the 9 pattern
    min_r, min_c = float('inf'), float('inf')
    max_r, max_c = -1, -1
    
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] == 9:
                min_r = min(min_r, i)
                max_r = max(max_r, i)
                min_c = min(min_c, j)
                max_c = max(max_c, j)
    
    # Extract the pattern
    pattern = []
    for i in range(min_r, max_r + 1):
        row = []
        for j in range(min_c, max_c + 1):
            row.append(grid[i][j])
        pattern.append(row)
    
    # Print the pattern
    for row in pattern:
        print(' '.join(map(str, row)))

# Parse input grid
input_str = """4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 9 4 4 9 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 9 4 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 9 9 9 9 4 4 9 4 4 4 4 4 4 4 4 4 4 4 4
4 4 9 9 4 9 9 9 9 4 4 4 3 4 4 4 4 4 4 4 4 4
4 4 4 4 4 9 4 4 9 4 4 4 3 4 3 4 4 4 4 4 4 4"""

grid = [[int(x) for x in row.split()] for row in input_str.strip().split('\n')]
find_9_pattern(grid)