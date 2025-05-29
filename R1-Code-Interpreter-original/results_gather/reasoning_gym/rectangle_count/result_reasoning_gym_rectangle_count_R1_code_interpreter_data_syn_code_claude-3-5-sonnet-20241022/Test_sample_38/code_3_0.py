def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_border(char):
        return char in '#█'
    
    def find_rectangle_from_corner(start_y, start_x):
        # Find width
        end_x = start_x
        while end_x < width and is_border(grid[start_y][end_x]):
            end_x += 1
        end_x -= 1
        
        # Find height
        end_y = start_y
        while end_y < height and is_border(grid[end_y][start_x]):
            end_y += 1
        end_y -= 1
        
        # Verify rectangle
        # Check all corners exist
        if not (is_border(grid[start_y][start_x]) and 
                is_border(grid[start_y][end_x]) and
                is_border(grid[end_y][start_x]) and
                is_border(grid[end_y][end_x])):
            return None
            
        # Check all edges are complete
        for x in range(start_x, end_x + 1):
            if not (is_border(grid[start_y][x]) and is_border(grid[end_y][x])):
                return None
        for y in range(start_y, end_y + 1):
            if not (is_border(grid[y][start_x]) and is_border(grid[y][end_x])):
                return None
                
        return (start_y, start_x, end_y, end_x)
    
    rectangles = set()
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Check if it's a potential top-left corner
                is_top = y == 0 or not is_border(grid[y-1][x])
                is_left = x == 0 or not is_border(grid[y][x-1])
                
                if is_top and is_left:
                    rect = find_rectangle_from_corner(y, x)
                    if rect:
                        rectangles.add(rect)
                        
                # Special case: check for overlapping rectangles at '█'
                if grid[y][x] == '█':
                    # Try finding rectangles in all directions from this point
                    for dy, dx in [(0,1), (1,0), (0,-1), (-1,0)]:
                        ny, nx = y + dy, x + dx
                        if 0 <= ny < height and 0 <= nx < width and is_border(grid[ny][nx]):
                            rect = find_rectangle_from_corner(min(y, ny), min(x, nx))
                            if rect:
                                rectangles.add(rect)
    
    return len(rectangles)

# Grid data (same as before)
grid = [
    # ... (same grid data as before)
]

print(find_rectangles(grid))