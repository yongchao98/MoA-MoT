def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def find_patterns(grid):
    sections = []
    current_section = []
    for row in grid:
        if all(x == 0 for x in row):
            if current_section:
                sections.append(current_section)
                current_section = []
        else:
            current_section.append(row)
    if current_section:
        sections.append(current_section)
    return sections

def get_column_pattern(section, col):
    # Extract non-zero values in the column
    values = [row[col] for row in section if row[col] != 0]
    if not values:
        return [0, 0, 0]
    
    # Analyze the pattern
    if len(values) == 1:
        return [values[0], 0, 0]
    elif len(values) == 2:
        # Check if it's part of a horizontal pattern
        for row in section:
            if sum(1 for x in row if x == values[0]) >= 3:
                return [values[0], values[0], values[0]]
        return [values[0], values[1], 0]
    else:
        # Check for special patterns
        if values.count(values[0]) >= 3:
            return [values[0], 0, values[0]]
        return [values[0], values[1], values[2]]

def transform_input(input_grid_str):
    grid = parse_grid(input_grid_str)
    sections = find_patterns(grid)
    width = len(grid[0])
    
    # Initialize output grid
    output = [[0]*width for _ in range(3)]
    
    # Process each section
    for section in sections:
        for col in range(width):
            pattern = get_column_pattern(section, col)
            # Update output if pattern contains non-zero values
            for row in range(3):
                if pattern[row] != 0:
                    output[row][col] = pattern[row]
    
    # Format and print result
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

transform_input(test_input)