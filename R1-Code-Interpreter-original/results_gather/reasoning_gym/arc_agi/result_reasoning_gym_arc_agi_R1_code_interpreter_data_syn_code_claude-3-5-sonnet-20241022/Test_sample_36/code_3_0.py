def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def get_vertical_pattern(section, col):
    # Get non-zero values in column
    values = [(i, row[col]) for i, row in enumerate(section) if row[col] != 0]
    if not values:
        return [0, 0, 0]
    
    num = values[0][1]
    
    # Pattern for number 9
    if num == 9:
        # Check for horizontal pattern of three 9s
        for i, row in enumerate(section):
            if row[col] == 9 and sum(1 for x in row if x == 9) >= 3:
                return [9, 0, 9]
        # Check for vertical pattern
        if len(values) >= 2:
            return [9, 9, 9]
        return [0, 9, 0]
    
    # Pattern for numbers 4, 5, 6
    # Check for horizontal groups of three
    for i, row in enumerate(section):
        if row[col] in [4, 5, 6]:
            count = sum(1 for x in row if x == row[col])
            if count >= 3:
                pos = len([x for x in row[:col] if x == row[col]])
                if pos == 0:
                    return [row[col], 0, row[col]]
                elif pos == 1:
                    return [0, row[col], 0]
                else:
                    return [row[col], 0, row[col]]
    
    # Default vertical pattern
    if len(values) >= 3:
        return [values[0][1], values[1][1], values[2][1]]
    elif len(values) == 2:
        return [values[0][1], values[1][1], 0]
    else:
        return [values[0][1], 0, 0]

def transform_grid(input_grid_str):
    grid = parse_grid(input_grid_str)
    width = len(grid[0])
    output = [[0]*width for _ in range(3)]
    
    # Process grid in sections
    current_section = []
    for row in grid:
        if all(x == 0 for x in row):
            if current_section:
                # Process current section
                for col in range(width):
                    pattern = get_vertical_pattern(current_section, col)
                    for i in range(3):
                        if pattern[i] != 0:
                            output[i][col] = pattern[i]
                current_section = []
        else:
            current_section.append(row)
    
    # Process last section if exists
    if current_section:
        for col in range(width):
            pattern = get_vertical_pattern(current_section, col)
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