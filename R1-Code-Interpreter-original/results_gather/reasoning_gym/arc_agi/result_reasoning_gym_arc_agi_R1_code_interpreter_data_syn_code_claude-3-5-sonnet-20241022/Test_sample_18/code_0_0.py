def find_patterns(grid):
    patterns = []
    rows = len(grid)
    cols = len(grid[0])
    
    # Convert input string to 2D grid
    grid = [[int(x) for x in row.split()] for row in grid]
    
    # Function to check if a 3x3 region contains a valid pattern
    def is_valid_pattern(r, c):
        if r + 2 >= rows or c + 2 >= cols:
            return False
        
        # Count 7s and their positions
        sevens = []
        for i in range(3):
            for j in range(3):
                if grid[r+i][c+j] == 7:
                    sevens.append((i, j))
        
        return len(sevens) >= 2  # Need at least 2 sevens to form a pattern
    
    # Find all patterns in the grid
    for i in range(rows-2):
        for j in range(cols-2):
            if is_valid_pattern(i, j):
                pattern = []
                for x in range(3):
                    row = []
                    for y in range(3):
                        if grid[i+x][j+y] == 7:
                            row.append(7)
                        else:
                            row.append(9)
                    pattern.append(row)
                patterns.append(pattern)
    
    # Generate output based on patterns
    output = []
    for pattern in patterns:
        # Convert pattern to output format
        out_pattern = []
        for row in pattern:
            if 7 in row:
                if row.count(7) >= 2:
                    out_pattern.append("7 7 7")
                else:
                    pos = row.index(7)
                    if pos == 0:
                        out_pattern.append("7 9 9")
                    elif pos == 1:
                        out_pattern.append("9 7 9")
                    else:
                        out_pattern.append("9 9 7")
            else:
                out_pattern.append("9 7 9")
        output.extend(out_pattern)
    
    return output

# Test input
test_input = [
    "9 9 9 9 2 9 9 9 9",
    "9 9 9 9 9 9 9 9 9",
    "9 9 9 9 9 9 9 9 9",
    "9 9 9 9 9 9 9 2 9",
    "9 9 9 9 9 9 9 9 9",
    "7 7 7 9 9 9 9 9 9",
    "7 9 7 9 9 9 9 9 9",
    "9 7 9 9 9 9 9 9 9",
    "9 9 9 9 9 9 9 9 9"
]

result = find_patterns(test_input)
for row in result:
    print(row)