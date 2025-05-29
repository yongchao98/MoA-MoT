def parse_grid(input_str):
    return [[int(x) for x in line.split()] for line in input_str.strip().split('\n')]

def find_pattern(input_grid):
    # Split the grid into three sections (divided by rows of 2's)
    sections = []
    current_section = []
    
    for row in input_grid:
        if all(x == 2 for x in row):
            if current_section:
                sections.append(current_section)
                current_section = []
        else:
            current_section.append(row)
    if current_section:
        sections.append(current_section)
    
    # Create 3x3 output grid
    output = [[3 for _ in range(3)] for _ in range(3)]
    
    # Process each section
    for section_idx, section in enumerate(sections):
        # For each column in output (0, 1, 2)
        for col_idx in range(3):
            # Count 5's in the corresponding region
            start_col = col_idx * 4  # 0, 4, 8
            
            # Count 5's in each row
            row_with_fives = set()
            for row_idx, row in enumerate(section):
                if any(row[col] == 5 for col in range(start_col, start_col + 3)):
                    row_with_fives.add(row_idx)
            
            # If we have exactly 2 rows containing 5's
            if len(row_with_fives) == 2:
                output[section_idx][col_idx] = 6
    
    return output

# Test input
test_input = """3 3 3 2 3 5 3 2 3 3 3
5 3 3 2 3 3 3 2 3 5 3
3 3 5 2 3 3 3 2 3 3 3
2 2 2 2 2 2 2 2 2 2 2
3 3 3 2 5 3 3 2 3 3 5
3 3 5 2 3 5 3 2 3 3 3
3 3 3 2 3 3 3 2 3 5 3
2 2 2 2 2 2 2 2 2 2 2
3 3 3 2 3 3 3 2 3 3 3
3 5 3 2 3 3 5 2 3 3 3
3 3 3 2 3 3 3 2 5 3 5"""

input_grid = parse_grid(test_input)
result = find_pattern(input_grid)

# Print result in required format
for row in result:
    print(' '.join(map(str, row)))