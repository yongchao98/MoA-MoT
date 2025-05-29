def find_sections(grid):
    # Convert input string to 2D grid
    rows = []
    valid_nums = {'0', '1', '3'}
    sections = []
    
    # Parse input grid
    for line in grid.strip().split('\n'):
        rows.append(line.split())
    
    # Find continuous sections of valid numbers
    def extract_section(start_row, start_col):
        section = []
        r, c = start_row, start_col
        while r < 10 and rows[r][c] != '9':
            row = []
            c = start_col
            while c < 10 and rows[r][c] != '9':
                if rows[r][c] in valid_nums:
                    row.append(rows[r][c])
                c += 1
            if row:  # Only add non-empty rows
                section.append(row)
            r += 1
        return section if section else None

    # Scan the grid for sections
    for i in range(10):
        for j in range(10):
            if rows[i][j] in valid_nums:
                # Check if this is the start of a new section
                if (i == 0 or rows[i-1][j] == '9') and (j == 0 or rows[i][j-1] == '9'):
                    section = extract_section(i, j)
                    if section:
                        sections.append(section)
    
    # Return the largest section (in terms of total cells)
    if sections:
        max_section = max(sections, key=lambda x: len(x) * len(x[0]))
        return max_section
    return None

# Test input
test_input = """9 9 9 9 9 9 9 9 9 9
9 9 9 9 9 3 0 1 3 9
9 9 9 9 9 0 3 0 0 9
9 0 0 3 9 1 0 1 3 9
9 3 0 1 9 9 9 9 9 9
9 1 0 0 9 9 9 9 9 9
9 0 3 1 9 9 0 0 0 0
9 0 0 0 9 9 3 0 1 3
9 3 0 1 9 9 0 3 0 0
9 9 9 9 9 9 0 0 0 1"""

result = find_sections(test_input)
if result:
    # Print the output grid
    for row in result:
        print(' '.join(row))