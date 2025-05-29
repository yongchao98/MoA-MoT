def find_sections(grid):
    # Convert input string to 2D grid
    rows = [line.split() for line in grid.strip().split('\n')]
    valid_nums = {'0', '1', '3'}
    sections = []
    
    def is_valid_start(r, c):
        # Check if position is valid start of a section
        if rows[r][c] not in valid_nums:
            return False
        # Check if there's a 9 above or to the left
        above = rows[r-1][c] if r > 0 else '9'
        left = rows[r][c-1] if c > 0 else '9'
        return above == '9' or left == '9'
    
    def get_section_size(r, c):
        # Get width and height of section
        width = height = 0
        i, j = r, c
        
        # Get height
        while i < 10 and rows[i][c] != '9':
            height += 1
            i += 1
            
        # Get width
        while j < 10 and rows[r][j] != '9':
            width += 1
            j += 1
            
        return height, width
    
    def extract_section(r, c, height, width):
        section = []
        for i in range(r, r + height):
            row = []
            for j in range(c, c + width):
                if rows[i][j] in valid_nums:
                    row.append(rows[i][j])
            if row:  # Only add non-empty rows
                section.append(row)
        return section
    
    # Find all valid sections
    for i in range(10):
        for j in range(10):
            if is_valid_start(i, j):
                height, width = get_section_size(i, j)
                section = extract_section(i, j, height, width)
                if section and all(len(row) == len(section[0]) for row in section):  # Ensure rectangular
                    sections.append(section)
    
    # Return the largest valid section
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
    for row in result:
        print(' '.join(row))