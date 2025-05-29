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
            
            # Count 5's in first and last rows of the section
            first_row_five = any(section[0][col] == 5 for col in range(start_col, start_col + 3))
            last_row_five = any(section[-1][col] == 5 for col in range(start_col, start_col + 3))
            
            # Count total 5's in this region
            total_fives = sum(1 for row in section for col in range(start_col, start_col + 3) if row[col] == 5)
            
            # If we have a 5 in both first and last row, or exactly 2 5's total
            if (first_row_five and last_row_five) or total_fives == 2:
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