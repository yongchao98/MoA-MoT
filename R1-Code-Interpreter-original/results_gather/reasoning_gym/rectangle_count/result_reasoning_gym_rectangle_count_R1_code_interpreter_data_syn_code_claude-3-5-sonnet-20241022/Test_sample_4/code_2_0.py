def count_rectangles(grid):
    # Quick count of obvious rectangles
    count = 0
    width = len(grid[0])
    height = len(grid)
    
    # Simple line scan
    for y in range(height):
        for x in range(width):
            # Only check positions that could be top-left corners
            if grid[y][x] in '#█':
                # Skip if not a potential corner
                if y > 0 and x > 0:
                    if grid[y-1][x] in '#█' or grid[y][x-1] in '#█':
                        continue
                
                # Quick check for minimum rectangle size
                if x + 1 >= width or y + 1 >= height:
                    continue
                if grid[y][x+1] not in '#█' or grid[y+1][x] not in '#█':
                    continue
                
                # Found a valid corner, count this rectangle
                count += 1
    
    print(count)

# Read input efficiently
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines at the end
            continue
        grid.append(line)
except EOFError:
    pass

count_rectangles(grid)