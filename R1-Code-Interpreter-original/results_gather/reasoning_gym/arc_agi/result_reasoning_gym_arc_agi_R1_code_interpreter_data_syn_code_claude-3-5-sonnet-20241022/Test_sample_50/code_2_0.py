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
        # Check left column (0-3)
        left_five = False
        for row in section:
            if row[0] == 5 or row[1] == 5 or row[2] == 5:
                left_five = True
                
        # Check middle column (4-7)
        middle_five = False
        for row in section:
            if row[4] == 5 or row[5] == 5 or row[6] == 5:
                middle_five = True
                
        # Check right column (8-10)
        right_five = False
        for row in section:
            if row[8] == 5 or row[9] == 5 or row[10] == 5:
                right_five = True
        
        # Update output based on 5's presence
        if left_five:
            output[section_idx][0] = 6
        if middle_five:
            output[section_idx][1] = 6
        if right_five:
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