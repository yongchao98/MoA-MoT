def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_border(char):
        return char in '#█'

    def find_complete_rectangle(top, left):
        # Find right boundary
        right = left
        while right < width and is_border(grid[top][right]):
            right += 1
        right -= 1
        
        # Find bottom boundary
        bottom = top
        while bottom < height and is_border(grid[bottom][left]):
            bottom += 1
        bottom -= 1
        
        # Validate rectangle
        if bottom <= top or right <= left:
            return None
            
        # Check if all borders are complete
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return None
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return None
                
        return (top, left, bottom, right)

    # Store found rectangles to avoid duplicates
    rectangles = set()
    
    # Scan for rectangle top-left corners
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Check if this could be a top-left corner
                if (x == 0 or not is_border(grid[y][x-1])) and \
                   (y == 0 or not is_border(grid[y-1][x])):
                    rect = find_complete_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)

    print(len(rectangles))

# Your ASCII grid
grid = """                                                                                
[... rest of the grid ...]
                                                                                """

find_rectangles(grid)