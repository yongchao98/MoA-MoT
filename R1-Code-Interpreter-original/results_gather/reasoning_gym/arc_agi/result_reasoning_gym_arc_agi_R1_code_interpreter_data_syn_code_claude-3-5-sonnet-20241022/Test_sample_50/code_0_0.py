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
    
    # For each section
    for section_idx, section in enumerate(sections):
        # Count 5's in each column
        for col in range(len(section[0])):
            fives_count = sum(1 for row in section if row[col] == 5)
            if fives_count >= 2:  # If there are 2 or more 5's in a column
                # Mark corresponding position in output as 6
                output[section_idx][col // 4] = 6
    
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