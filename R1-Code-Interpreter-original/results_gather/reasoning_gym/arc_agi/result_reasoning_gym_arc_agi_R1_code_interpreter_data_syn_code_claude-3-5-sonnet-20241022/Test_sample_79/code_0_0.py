def parse_grid(input_grid):
    # Remove the border rows and columns that contain 7s
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the actual content boundaries (excluding 7s)
    start_row = 0
    end_row = rows
    start_col = 0
    end_col = cols
    
    for i in range(rows):
        if any(x != 7 for x in input_grid[i]):
            start_row = i
            break
    
    for i in range(rows-1, -1, -1):
        if any(x != 7 for x in input_grid[i]):
            end_row = i + 1
            break
    
    for j in range(cols):
        if any(input_grid[i][j] != 7 for i in range(rows)):
            start_col = j
            break
    
    for j in range(cols-1, -1, -1):
        if any(input_grid[i][j] != 7 for i in range(rows)):
            end_col = j + 1
            break
    
    # Find horizontal sections
    sections = []
    current_row = start_row
    while current_row < end_row:
        # Find the height of this section
        section_height = 1
        while (current_row + section_height < end_row and 
               all(input_grid[current_row][j] == input_grid[current_row + section_height][j] 
                   for j in range(start_col, end_col))):
            section_height += 1
        
        # Get the distinct numbers in this section from left to right
        row = input_grid[current_row][start_col:end_col]
        distinct_nums = []
        prev = None
        for num in row:
            if num != prev and num != 7:
                distinct_nums.append(num)
                prev = num
        
        sections.append(distinct_nums)
        current_row += section_height
    
    # Format the output
    result = ""
    for section in sections:
        result += " ".join(map(str, section)) + "\n"
    
    print(result.strip())

# Test input
test_input = [
    [7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7],
    [7,2,2,2,2,2,2,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,2,2,2,2,2,2,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,2,2,2,2,2,2,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,2,2,2,2,2,2,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,2,2,2,2,2,2,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,5,5,5,5,5,5,9,9,9,9,9,9,9,9,5,5,5,5,7,7],
    [7,0,0,0,0,0,0,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,0,0,0,0,0,0,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,0,0,0,0,0,0,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,0,0,0,0,0,0,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,0,0,0,0,0,0,6,6,6,6,6,6,6,6,2,2,2,2,7,7],
    [7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]
]

parse_grid(test_input)