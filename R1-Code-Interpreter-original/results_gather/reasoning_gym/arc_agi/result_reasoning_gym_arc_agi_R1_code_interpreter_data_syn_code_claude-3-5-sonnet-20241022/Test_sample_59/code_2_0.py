def transform_grid(input_grid):
    n = len(input_grid)
    output = [[0]*n for _ in range(n)]
    
    def find_pattern_bounds():
        patterns = []
        for i in range(n):
            for j in range(n):
                if input_grid[i][j] != 0:
                    # Check if this is the top-left corner of a pattern
                    if (i == 0 or input_grid[i-1][j] == 0) and (j == 0 or input_grid[i][j-1] == 0):
                        # Find pattern bounds
                        height = 0
                        width = 0
                        while i + height < n and input_grid[i + height][j] != 0:
                            height += 1
                        while j + width < n and input_grid[i][j + width] != 0:
                            width += 1
                        patterns.append((i, j, height, width))
        return patterns
    
    def expand_pattern(top, left, height, width):
        frame_num = input_grid[top][left]
        inner_num = input_grid[top+1][left+1] if height > 2 and width > 2 else frame_num
        
        # New bounds
        new_top = max(0, top - 1)
        new_left = max(0, left - 1)
        new_bottom = min(n-1, top + height)
        new_right = min(n-1, left + width)
        
        # Fill outer frame
        for i in range(new_top, new_bottom + 1):
            for j in range(new_left, new_right + 1):
                if i == new_top or i == new_bottom or j == new_left or j == new_right:
                    output[i][j] = frame_num
                elif i >= top and i < top + height and j >= left and j < left + width:
                    # Preserve original pattern structure
                    output[i][j] = input_grid[i][j]
                else:
                    # Fill expanded inner area
                    output[i][j] = inner_num
    
    # Find and transform each pattern
    patterns = find_pattern_bounds()
    for top, left, height, width in patterns:
        expand_pattern(top, left, height, width)
    
    return output

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,6,4,6,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0]
]

result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))