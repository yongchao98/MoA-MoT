def find_patterns(grid):
    patterns = {}
    rows = len(grid)
    cols = len(grid[0])
    
    # Function to extract 3x3 pattern centered at (i,j)
    def get_pattern(i, j):
        pattern = []
        for di in [-1, 0, 1]:
            row = []
            for dj in [-1, 0, 1]:
                if 0 <= i+di < rows and 0 <= j+dj < cols:
                    row.append(grid[i+di][j+dj])
                else:
                    row.append(7)  # Default value
            pattern.append(row)
        return tuple(map(tuple, pattern))
    
    # Find all non-7 values and their surrounding patterns
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 7:
                pattern = get_pattern(i, j)
                if pattern not in patterns:
                    patterns[pattern] = 0
                patterns[pattern] += 1
    
    # Return the most common pattern
    if patterns:
        return max(patterns.items(), key=lambda x: x[1])[0]
    return None

# Parse the input grid
input_str = """7 7 7 7 7 7 7 7 7 9 9 7 7 7
7 7 7 7 7 7 7 7 9 9 9 7 7 7
7 7 7 7 7 7 7 7 7 9 7 7 7 7
7 7 5 7 5 7 7 7 7 7 7 7 7 7
7 7 7 5 7 7 7 7 7 7 7 7 9 9
7 7 5 5 5 7 7 7 7 7 7 9 9 9
7 7 7 7 7 7 7 7 8 7 7 7 9 7
7 7 7 7 7 7 7 8 7 8 7 7 7 7
7 7 7 7 7 7 7 7 8 7 7 7 7 7
7 7 7 9 9 7 7 7 7 7 7 5 7 5
7 7 9 9 9 7 7 7 7 7 7 7 5 7
7 7 7 9 7 7 7 7 7 7 7 5 5 5
7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 7 7 7 7 7 7 7 7 7 7 7"""

grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]

# Find the most common pattern
pattern = find_patterns(grid)

# Print the result
for row in pattern:
    print(' '.join(map(str, row)))