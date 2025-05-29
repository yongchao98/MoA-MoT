def analyze_grid(grid):
    # Convert input string to 2D grid
    grid = [[int(x) for x in row.split()] for row in grid]
    rows, cols = len(grid), len(grid[0])
    
    # Find connected groups of 7s
    def get_pattern(r, c):
        pattern = []
        # Check 3x3 area
        for i in range(3):
            if r+i < rows:
                row = []
                for j in range(3):
                    if c+j < cols:
                        row.append(grid[r+i][c+j])
                    else:
                        row.append(9)
                pattern.append(row)
            else:
                pattern.append([9, 9, 9])
        return pattern

    # Find significant patterns of 7s
    patterns = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 7:
                pattern = get_pattern(i, j)
                # Count 7s in pattern
                seven_count = sum(row.count(7) for row in pattern)
                if seven_count >= 2:  # Only consider patterns with 2 or more 7s
                    patterns.append(pattern)

    # If no patterns found, return empty result
    if not patterns:
        return []

    # Get the most significant pattern (the one with most 7s)
    pattern = patterns[0]
    
    # Generate output based on pattern
    output = []
    for row in pattern:
        seven_count = row.count(7)
        if seven_count == 0:
            output.append("9 7 9")
        elif seven_count == 1:
            if row[0] == 7:
                output.append("7 9 9")
            elif row[1] == 7:
                output.append("9 7 9")
            else:
                output.append("9 9 7")
        elif seven_count >= 2:
            output.append("7 7 7")
    
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

result = analyze_grid(test_input)
for row in result:
    print(row)