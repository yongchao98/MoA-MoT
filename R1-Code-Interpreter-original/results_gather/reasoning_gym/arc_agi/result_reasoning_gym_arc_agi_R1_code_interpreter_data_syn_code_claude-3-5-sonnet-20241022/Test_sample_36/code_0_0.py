def parse_grid(grid_str):
    return [[int(x) for x in row.split()] for row in grid_str.strip().split('\n')]

def find_sections(grid):
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
    return sections

def transform_section(section):
    # Analyze patterns in each section
    result = []
    height = len(section)
    width = len(section[0])
    
    # Create three rows for output
    for i in range(3):
        new_row = [0] * width
        for col in range(width):
            # Look for non-zero numbers in this column
            numbers = [row[col] for row in section if row[col] != 0]
            if numbers:
                if i == 0 and len(numbers) >= 1:
                    new_row[col] = numbers[0]
                elif i == 1 and len(numbers) >= 2:
                    new_row[col] = numbers[1]
                elif i == 2 and len(numbers) >= 3:
                    new_row[col] = numbers[2]
        result.append(new_row)
    return result

def predict_output(input_grid_str):
    # Parse input grid
    grid = parse_grid(input_grid_str)
    
    # Find sections between zero rows
    sections = find_sections(grid)
    
    # Transform each section and combine results
    final_output = []
    for section in sections:
        transformed = transform_section(section)
        if not final_output:
            final_output = transformed
        else:
            # Combine patterns from multiple sections
            for i in range(3):
                for j in range(len(transformed[0])):
                    if transformed[i][j] != 0:
                        final_output[i][j] = transformed[i][j]
    
    # Format output
    result = ""
    for row in final_output:
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

predict_output(test_input)