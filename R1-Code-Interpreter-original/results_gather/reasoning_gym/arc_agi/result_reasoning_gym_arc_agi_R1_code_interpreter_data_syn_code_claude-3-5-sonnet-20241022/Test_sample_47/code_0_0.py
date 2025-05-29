def parse_grid(input_str):
    return [[int(x) for x in line.split()] for line in input_str.strip().split('\n')]

def count_ones_in_sections(grid):
    # Split grid into three sections and count 1s in each
    sections = []
    current_section = []
    section_counts = []
    
    for row in grid:
        if 9 in row:  # This is a divider row
            if current_section:
                sections.append(current_section)
                section_counts.append(sum(x == 1 for row in current_section for x in row))
                current_section = []
        else:
            current_section.append(row)
    
    # Add the last section
    if current_section:
        sections.append(current_section)
        section_counts.append(sum(x == 1 for row in current_section for x in row))
    
    return section_counts

def analyze_pattern(section_counts):
    # Convert counts to output numbers based on pattern
    # If there are 2 or more 1s in a section, output 4
    # Otherwise output 6
    return [4 if count >= 2 else 6 for count in section_counts]

def format_output(numbers):
    return '\n'.join(' '.join(str(x) for x in row) for row in [numbers])

# Test input
test_input = """6 6 6 9 6 6 6 9 6 1 6
6 1 6 9 6 6 6 9 6 6 6
6 6 6 9 6 1 6 9 1 6 6
9 9 9 9 9 9 9 9 9 9 9
6 6 6 9 6 6 1 9 6 6 6
6 6 6 9 6 1 6 9 6 6 1
6 1 6 9 6 6 6 9 6 6 6
9 9 9 9 9 9 9 9 9 9 9
1 6 6 9 6 6 6 9 6 6 6
6 6 6 9 1 6 6 9 6 1 6
1 6 6 9 6 6 1 9 6 6 6"""

grid = parse_grid(test_input)
section_counts = count_ones_in_sections(grid)
output_numbers = analyze_pattern(section_counts)
print(format_output(output_numbers))