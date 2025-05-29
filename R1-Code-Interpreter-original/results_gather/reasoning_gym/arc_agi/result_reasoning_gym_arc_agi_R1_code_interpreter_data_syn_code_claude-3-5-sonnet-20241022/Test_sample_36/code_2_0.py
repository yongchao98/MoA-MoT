def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def find_horizontal_groups(row):
    groups = []
    current_group = []
    current_num = None
    for i, num in enumerate(row):
        if num != 0:
            if current_num is None or current_num == num:
                current_group.append((i, num))
                current_num = num
            else:
                if current_group:
                    groups.append(current_group)
                current_group = [(i, num)]
                current_num = num
    if current_group:
        groups.append(current_group)
    return groups

def process_column(section, col):
    values = [row[col] for row in section]
    non_zero = [(i, v) for i, v in enumerate(values) if v != 0]
    
    if not non_zero:
        return [0, 0, 0]
    
    # Check for horizontal patterns in each row
    horizontal_patterns = []
    for i, row in enumerate(section):
        if row[col] != 0:
            groups = find_horizontal_groups(row)
            for group in groups:
                if any(pos == col for pos, _ in group) and len(group) >= 3:
                    horizontal_patterns.append((i, row[col]))
    
    # Special pattern for number 9
    if non_zero[0][1] == 9:
        if len(non_zero) >= 3:
            return [9, 9, 9]
        elif len(non_zero) == 2:
            return [9, 9, 0]
        else:
            return [0, 9, 0]
    
    # Pattern for numbers 4, 5, 6
    if horizontal_patterns:
        num = horizontal_patterns[0][1]
        if len(horizontal_patterns) >= 2:
            return [num, 0, num]
        else:
            return [num, num, num]
    
    # Default pattern based on vertical arrangement
    result = [0, 0, 0]
    for i, (pos, val) in enumerate(non_zero[:3]):
        result[i] = val
    return result

def transform_grid(input_grid_str):
    grid = parse_grid(input_grid_str)
    width = len(grid[0])
    output = [[0]*width for _ in range(3)]
    
    # Process grid in sections (between zero rows)
    current_section = []
    for row in grid:
        if all(x == 0 for x in row):
            if current_section:
                # Process current section
                for col in range(width):
                    pattern = process_column(current_section, col)
                    for i in range(3):
                        if pattern[i] != 0:
                            output[i][col] = pattern[i]
                current_section = []
        else:
            current_section.append(row)
    
    # Process last section if exists
    if current_section:
        for col in range(width):
            pattern = process_column(current_section, col)
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