def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def get_horizontal_group(row, col):
    # Get the size and position of horizontal group containing this position
    num = row[col]
    if num == 0:
        return 0, 0
    
    # Count consecutive same numbers left and right
    left = 0
    for i in range(col-1, -1, -1):
        if row[i] == num:
            left += 1
        else:
            break
            
    right = 0
    for i in range(col+1, len(row)):
        if row[i] == num:
            right += 1
        else:
            break
            
    return left + right + 1, left

def process_column(section, col):
    values = [(i, row[col]) for i, row in enumerate(section) if row[col] != 0]
    if not values:
        return [0, 0, 0]
        
    # Process number 9
    if values[0][1] == 9:
        # Check for horizontal pattern
        for i, row in enumerate(section):
            group_size, pos = get_horizontal_group(row, col)
            if group_size == 3:
                if pos == 0:
                    return [9, 0, 9]
                elif pos == 1:
                    return [0, 9, 0]
                else:
                    return [9, 0, 9]
        # Vertical pattern
        if len(values) >= 2:
            return [9, 9, 9]
        return [0, 9, 0]
    
    # Process numbers 4, 5, 6
    num = values[0][1]
    # Check for horizontal groups
    for i, row in enumerate(section):
        if row[col] in [4, 5, 6]:
            group_size, pos = get_horizontal_group(row, col)
            if group_size == 3:
                if pos == 0:
                    return [num, 0, num]
                elif pos == 1:
                    return [0, num, 0]
                else:
                    return [num, 0, num]
    
    # Default pattern based on vertical arrangement
    if len(values) == 1:
        return [num, 0, 0]
    elif len(values) == 2:
        return [values[0][1], values[1][1], 0]
    else:
        return [values[0][1], values[1][1], values[2][1]]

def transform_grid(input_grid_str):
    grid = parse_grid(input_grid_str)
    width = len(grid[0])
    output = [[0]*width for _ in range(3)]
    
    # Find sections between zero rows
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
    
    # Process each section
    for section in sections:
        for col in range(width):
            pattern = process_column(section, col)
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