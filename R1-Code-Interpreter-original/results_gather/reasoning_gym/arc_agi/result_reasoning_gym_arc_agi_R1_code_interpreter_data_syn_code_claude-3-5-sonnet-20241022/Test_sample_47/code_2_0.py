def parse_grid(input_str):
    return [[int(x) for x in line.split()] for line in input_str.strip().split('\n')]

def analyze_section(section):
    # Split section into three parts (left, middle, right) based on 9s
    left = [row[:3] for row in section]
    middle = [row[4:7] for row in section]
    right = [row[8:] for row in section]
    
    # Count 1s in each part
    left_ones = sum(x == 1 for row in left for x in row)
    middle_ones = sum(x == 1 for row in middle for x in row)
    right_ones = sum(x == 1 for row in right for x in row)
    
    # Generate output row based on number of 1s in each part
    return [4 if left_ones >= 2 else 6, 
            4 if middle_ones >= 2 else 6, 
            4 if right_ones >= 2 else 6]

def process_grid(grid):
    # Split into three main sections
    section1 = grid[0:3]
    section2 = grid[4:7]
    section3 = grid[8:11]
    
    # Analyze each section
    output = []
    for section in [section1, section2, section3]:
        output.append(analyze_section(section))
    return output

def format_output(output_grid):
    return '\n'.join(' '.join(str(x) for x in row) for row in output_grid)

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
output_grid = process_grid(grid)
print(format_output(output_grid))