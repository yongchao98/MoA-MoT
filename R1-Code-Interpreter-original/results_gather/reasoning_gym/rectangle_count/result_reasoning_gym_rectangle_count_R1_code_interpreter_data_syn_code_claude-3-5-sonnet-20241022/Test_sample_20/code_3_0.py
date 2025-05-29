def quick_count_rectangles():
    # Process input more efficiently
    grid = []
    chars = set(['#', '█'])
    first_char_col = float('inf')
    last_char_col = -1
    first_char_row = float('inf')
    last_char_row = -1
    
    # Read and process input, tracking boundaries
    row = 0
    while True:
        try:
            line = input()
            if any(c in chars for c in line):
                grid.append(line)
                # Track boundaries
                for col, c in enumerate(line):
                    if c in chars:
                        first_char_col = min(first_char_col, col)
                        last_char_col = max(last_char_col, col)
                        first_char_row = min(first_char_row, row)
                        last_char_row = max(last_char_row, row)
                row += 1
        except EOFError:
            break
    
    if not grid:
        return 0
    
    # Trim grid to relevant area
    grid = [row[first_char_col:last_char_col+1] for row in grid[first_char_row:last_char_row+1]]
    
    count = 0
    height = len(grid)
    width = len(grid[0]) if grid else 0
    
    # Fast rectangle detection
    for i in range(height-1):
        for j in range(width-1):
            if grid[i][j] in chars:
                # Check if this could be top-left corner
                if (i == 0 or grid[i-1][j] not in chars) and (j == 0 or grid[i][j-1] not in chars):
                    # Find right edge
                    right = j + 1
                    while right < width and grid[i][right] in chars:
                        right += 1
                    right -= 1
                    
                    if right > j:  # Found potential width
                        # Find bottom edge
                        bottom = i + 1
                        while bottom < height and grid[bottom][j] in chars:
                            bottom += 1
                        bottom -= 1
                        
                        if bottom > i:  # Found potential height
                            # Quick validate rectangle
                            valid = True
                            for x in range(i, bottom+1):
                                if grid[x][j] not in chars or grid[x][right] not in chars:
                                    valid = False
                                    break
                            if valid:
                                for y in range(j, right+1):
                                    if grid[i][y] not in chars or grid[bottom][y] not in chars:
                                        valid = False
                                        break
                            if valid:
                                count += 1
    
    return count

print(quick_count_rectangles())