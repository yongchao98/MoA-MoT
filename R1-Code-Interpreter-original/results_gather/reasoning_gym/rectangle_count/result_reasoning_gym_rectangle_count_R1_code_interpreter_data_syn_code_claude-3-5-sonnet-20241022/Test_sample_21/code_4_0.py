def count_rectangles(grid):
    # Count basic rectangles and overlaps in a single pass
    basic_rects = 0
    overlaps = 0
    
    # Process each character only once
    for y in range(len(grid)):
        for x in range(len(grid[y])):
            if grid[y][x] == '#':
                # Check if it's a potential top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    # Quick check if there's a rectangle
                    if x+1 < len(grid[y]) and y+1 < len(grid) and grid[y][x+1] in '#█' and grid[y+1][x] in '#█':
                        basic_rects += 1
            elif grid[y][x] == '█':
                overlaps += 1
    
    # Print total (basic rectangles + overlaps)
    print(basic_rects + overlaps)

# Read input
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

count_rectangles(grid)