def find_bottom_right_section(grid):
    # Convert input string to 2D grid
    rows = [line.split() for line in grid.strip().split('\n')]
    valid_nums = {'0', '1', '3'}
    
    # Find the middle point of the grid
    mid_row = len(rows) // 2
    mid_col = len(rows[0]) // 2
    
    def extract_section(r, c):
        section = []
        current_row = []
        
        # Extract numbers until we hit a 9 or grid boundary
        while r < len(rows) and rows[r][c] != '9':
            current_row = []
            j = c
            while j < len(rows[0]) and rows[r][j] != '9':
                if rows[r][j] in valid_nums:
                    current_row.append(rows[r][j])
                j += 1
            if current_row:  # Only add non-empty rows
                section.append(current_row)
            r += 1
        return section if section else None
    
    # Scan the bottom-right quadrant
    for i in range(mid_row, len(rows)):
        for j in range(mid_col, len(rows[0])):
            if rows[i][j] in valid_nums:
                # Check if this is the start of a section
                if (i == 0 or rows[i-1][j] == '9') and (j == 0 or rows[i][j-1] == '9'):
                    section = extract_section(i, j)
                    if section:
                        return section
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

result = find_bottom_right_section(test_input)
if result:
    for row in result:
        print(' '.join(row))