def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [row[:] for row in input_grid]
    
    # Find the middle section of zeros
    zero_section_start = None
    for i in range(rows):
        if all(x == 0 for x in input_grid[i]):
            if zero_section_start is None:
                zero_section_start = i
            zero_section_end = i
    
    # Process each column
    for col in range(cols):
        # Process upper half - find topmost special number
        top_special = -1
        for row in range(zero_section_start):
            if input_grid[row][col] not in [0, 5]:
                if top_special == -1:
                    top_special = row
        if top_special != -1:
            for r in range(top_special, zero_section_start):
                output[r][col] = 0
        
        # Process lower half - find bottommost special number
        bottom_special = -1
        for row in range(rows-1, zero_section_end, -1):
            if input_grid[row][col] not in [0, 5]:
                if bottom_special == -1:
                    bottom_special = row
        if bottom_special != -1:
            for r in range(zero_section_end + 1, bottom_special + 1):
                output[r][col] = 0
    
    return output

# Test input
test_input = [
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,6],
    [5,5,5,5,5,5,6,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,6,5,5],
    [5,6,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,6,5,5,5,5,5,5,5],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,6,5,5,5,5,5,5,5,6,5,5],
    [5,6,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,6,5,5,5,5,5,5,6,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5]
]

result = process_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))