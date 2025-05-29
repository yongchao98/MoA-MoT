def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Found top-left corner
                if (x == 0 or not is_border(grid[y][x-1])) and (y == 0 or not is_border(grid[y-1][x])):
                    # Find width
                    w = 1
                    while x + w < width and is_border(grid[y][x+w]):
                        w += 1
                    
                    # Find height
                    h = 1
                    while y + h < height and is_border(grid[y+h][x]):
                        h += 1
                    
                    # Verify rectangle
                    is_valid = True
                    for cy in range(y, y+h):
                        if not is_border(grid[cy][x]) or not is_border(grid[cy][x+w-1]):
                            is_valid = False
                            break
                    for cx in range(x, x+w):
                        if not is_border(grid[y][cx]) or not is_border(grid[y+h-1][cx]):
                            is_valid = False
                            break
                    
                    if is_valid:
                        rectangles.append((x, y, w, h))
    
    print(len(rectangles))

# Your grid
grid = """                                                                                
[... rest of the grid ...]
                                                                                """

find_rectangles(grid)