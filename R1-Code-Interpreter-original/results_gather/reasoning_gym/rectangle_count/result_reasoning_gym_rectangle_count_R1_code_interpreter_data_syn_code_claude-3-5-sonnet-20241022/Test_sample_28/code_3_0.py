def count_rectangles(grid):
    # Remove empty lines and get actual dimensions
    grid = [line for line in grid if line.strip()]
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    
    # Count base rectangles by scanning horizontal lines
    rectangles = 0
    overlaps = 0
    
    # Scan each line for start of horizontal edges
    for y in range(height):
        x = 0
        while x < width:
            # Skip non-border characters
            if grid[y][x] not in '#█':
                x += 1
                continue
            
            # Found start of horizontal line
            if x == 0 or grid[y][x-1] not in '#█':
                # Find end of horizontal line
                start_x = x
                while x < width and grid[y][x] in '#█':
                    x += 1
                end_x = x - 1
                
                # Check if this is top edge of a rectangle
                if y == 0 or grid[y-1][start_x] not in '#█':
                    # Find bottom edge
                    bottom_y = y + 1
                    is_rectangle = False
                    while bottom_y < height:
                        # Check if this is a complete bottom edge
                        if all(grid[bottom_y][i] in '#█' for i in range(start_x, end_x + 1)):
                            # Verify vertical edges
                            if all(grid[i][start_x] in '#█' and grid[i][end_x] in '#█' 
                                  for i in range(y, bottom_y + 1)):
                                rectangles += 1
                                is_rectangle = True
                        bottom_y += 1
            else:
                x += 1
    
    # Count overlaps (█ characters that are part of multiple rectangles)
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                # Check if this overlap point connects different rectangles
                if (x > 0 and x < width-1 and y > 0 and y < height-1 and
                    ((grid[y-1][x] in '#█' and grid[y+1][x] in '#█' and
                      grid[y][x-1] in '#█' and grid[y][x+1] in '#█'))):
                    overlaps += 1
    
    return rectangles + overlaps

# Read input
grid = []
try:
    while True:
        line = input()
        grid.append(line)
except EOFError:
    pass

print(count_rectangles(grid))