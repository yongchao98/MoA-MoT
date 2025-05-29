def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def find_sections(grid):
    sections = []
    current = []
    for row in grid:
        if all(x == 0 for x in row):
            if current:
                sections.append(current)
                current = []
        else:
            current.append(row)
    if current:
        sections.append(current)
    return sections

def get_pattern(section, col):
    # Get non-zero values in column
    values = [(i, row[col]) for i, row in enumerate(section) if row[col] != 0]
    if not values:
        return [0, 0, 0]
    
    # Check horizontal patterns for each row
    for i, row in enumerate(section):
        if row[col] != 0:
            # Count consecutive same numbers
            num = row[col]
            count = 0
            start_pos = 0
            for j, val in enumerate(row):
                if val == num:
                    if count == 0:
                        start_pos = j
                    count += 1
                elif count > 0:
                    break
            
            if count >= 3:  # Found horizontal group of 3 or more
                if col == start_pos:  # First position
                    return [num, 0, num] if num != 9 else [num, num, num]
                elif col == start_pos + 1:  # Middle position
                    return [0, num, 0] if num != 9 else [num, 0, num]
                else:  # Last position
                    return [num, 0, num] if num != 9 else [num, num, num]
    
    # Handle vertical patterns
    if values[0][1] == 9:
        if len(values) >= 2:
            return [9, 9, 9]
        return [0, 9, 0]
    
    # Default pattern
    result = [0, 0, 0]
    for i, (_, val) in enumerate(values[:3]):
        result[i] = val
    return result

def transform_grid(input_grid_str):
    grid = parse_grid(input_grid_str)
    sections = find_sections(grid)
    width = len(grid[0])
    output = [[0]*width for _ in range(3)]
    
    for section in sections:
        for col in range(width):
            pattern = get_pattern(section, col)
            for i in range(3):
                if pattern[i] != 0:
                    output[i][col] = pattern[i]
    
    # Print result
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = """0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
9 9 0 0 4 0 4 0 0 5 0 0 6 6 6
0 0 9 0 4 0 4 0 5 5 5 0 0 6 0
9 9 0 0 4 0 4 0 0 5 0 0 6 6 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
9 9 9 0 4 4 0 0 5 0 5 0 0 6 0
0 9 0 0 0 0 4 0 5 0 5 0 6 6 6
9 9 0 0 4 4 0 0 5 0 5 0 0 6 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"""

transform_grid(test_input)