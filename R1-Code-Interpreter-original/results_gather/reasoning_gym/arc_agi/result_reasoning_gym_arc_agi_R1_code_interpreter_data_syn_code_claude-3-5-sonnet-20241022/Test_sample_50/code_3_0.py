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
        # Count 5's in specific positions
        # Left region (columns 0-2)
        left_count = sum(1 for row in section for col in [0,1,2] if row[col] == 5)
        
        # Middle region (columns 4-6)
        middle_count = sum(1 for row in section for col in [4,5,6] if row[col] == 5)
        
        # Right region (columns 8-10)
        right_count = sum(1 for row in section for col in [8,9,10] if row[col] == 5)
        
        # Update output based on count of 5's
        # Only mark 6 if there are exactly 2 or more 5's in the region
        if left_count >= 2:
            output[section_idx][0] = 6
        if middle_count >= 2:
            output[section_idx][1] = 6
        if right_count >= 2:
            output[section_idx][2] = 6
    
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