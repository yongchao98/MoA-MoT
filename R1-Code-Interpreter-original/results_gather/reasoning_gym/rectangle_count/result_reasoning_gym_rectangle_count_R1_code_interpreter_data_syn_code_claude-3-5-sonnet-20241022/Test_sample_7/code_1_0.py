def find_rectangles(grid):
    # Convert grid to list of strings, removing empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_border(x, y):
        return 0 <= y < height and 0 <= x < width and grid[y][x] in '#█'
    
    def validate_rectangle(top, left, bottom, right):
        # Check dimensions
        if right <= left or bottom <= top:
            return False
            
        # Check if this is a minimal rectangle
        has_interior = False
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] == ' ':
                    has_interior = True
                    break
            if has_interior:
                break
        if not has_interior:
            return False
            
        # Check all borders
        # Top and bottom edges
        for x in range(left, right + 1):
            if not (is_border(x, top) and is_border(x, bottom)):
                return False
        
        # Left and right edges
        for y in range(top, bottom + 1):
            if not (is_border(left, y) and is_border(right, y)):
                return False
                
        return True

    rectangles = set()
    # Find potential corners (marked with # or █)
    corners = [(x, y) for y in range(height) for x in range(width) 
              if grid[y][x] in '#█']
    
    # Try all possible corner combinations
    for top_left in corners:
        for bottom_right in corners:
            if bottom_right[1] <= top_left[1] or bottom_right[0] <= top_left[0]:
                continue
            if validate_rectangle(top_left[1], top_left[0], 
                                bottom_right[1], bottom_right[0]):
                rectangles.add((top_left[0], top_left[1], 
                              bottom_right[0], bottom_right[1]))

    return len(rectangles)

grid = """                                                                                
[... your grid here ...]
                                                                                """

print(find_rectangles(grid))