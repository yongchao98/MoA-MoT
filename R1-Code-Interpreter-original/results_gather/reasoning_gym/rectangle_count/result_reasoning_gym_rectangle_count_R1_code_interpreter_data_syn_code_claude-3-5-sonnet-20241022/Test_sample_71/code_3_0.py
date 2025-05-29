def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_border(ch):
        return ch in '#█'
    
    def find_complete_rectangle(start_x, start_y):
        if not is_border(grid[start_y][start_x]):
            return None
            
        # Find right edge
        right = start_x
        while right < width and is_border(grid[start_y][right]):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = start_y
        while bottom < height and is_border(grid[bottom][start_x]):
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        if right <= start_x or bottom <= start_y:
            return None
            
        # Check all borders are complete
        for x in range(start_x, right + 1):
            if not (is_border(grid[start_y][x]) and is_border(grid[bottom][x])):
                return None
                
        for y in range(start_y, bottom + 1):
            if not (is_border(grid[y][start_x]) and is_border(grid[y][right])):
                return None
        
        # Check if it's not just a line
        has_space = False
        for y in range(start_y + 1, bottom):
            for x in range(start_x + 1, right):
                if grid[y][x] == ' ':
                    has_space = True
                    break
            if has_space:
                break
        
        if not has_space:
            return None
            
        return (start_x, start_y, right, bottom)
    
    rectangles = set()
    processed = set()
    
    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            if (x, y) in processed:
                continue
                
            rect = find_complete_rectangle(x, y)
            if rect:
                x1, y1, x2, y2 = rect
                # Mark all border points as processed
                for i in range(x1, x2 + 1):
                    processed.add((i, y1))
                    processed.add((i, y2))
                for j in range(y1, y2 + 1):
                    processed.add((x1, j))
                    processed.add((x2, j))
                rectangles.add(rect)
    
    # Filter out rectangles that are just parts of larger ones
    final_rectangles = set()
    for rect in rectangles:
        x1, y1, x2, y2 = rect
        is_unique = True
        for other in rectangles:
            if rect == other:
                continue
            ox1, oy1, ox2, oy2 = other
            # Check if this rectangle is completely contained within another
            if (ox1 <= x1 and x2 <= ox2 and oy1 <= y1 and y2 <= oy2):
                # Check if it's not just sharing a border
                if (ox1 < x1 or x2 < ox2) and (oy1 < y1 or y2 < oy2):
                    is_unique = False
                    break
        if is_unique:
            final_rectangles.add(rect)
    
    return len(final_rectangles)

# Test grid (your provided grid here)
test_grid = """                                                                                
[... rest of your grid ...]
             ###################                                                """

print(find_rectangles(test_grid))