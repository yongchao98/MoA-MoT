def find_rectangles(grid):
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    rectangles = set()

    def is_valid_char(c):
        return c in '#█'

    def trace_rectangle(start_y, start_x):
        # Find right edge
        right_x = start_x
        while right_x < width and is_valid_char(grid[start_y][right_x]):
            right_x += 1
        right_x -= 1

        # Find bottom edge
        bottom_y = start_y
        while bottom_y < height and is_valid_char(grid[bottom_y][start_x]):
            bottom_y += 1
        bottom_y -= 1

        # Verify rectangle
        for y in range(start_y, bottom_y + 1):
            if not is_valid_char(grid[y][start_x]) or not is_valid_char(grid[y][right_x]):
                return None
        for x in range(start_x, right_x + 1):
            if not is_valid_char(grid[start_y][x]) or not is_valid_char(grid[bottom_y][x]):
                return None
        
        return (start_y, start_x, bottom_y, right_x)

    # Find all potential rectangles
    for y in range(height):
        for x in range(width):
            if is_valid_char(grid[y][x]):
                # Check if this could be the top-left corner of a rectangle
                if (y == 0 or not is_valid_char(grid[y-1][x])) and \
                   (x == 0 or not is_valid_char(grid[y][x-1])):
                    rect = trace_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)

    # Count overlapping rectangles
    overlaps = set()
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                # Find all rectangles containing this point
                containing_rects = []
                for rect in rectangles:
                    if rect[0] <= y <= rect[2] and rect[1] <= x <= rect[3]:
                        containing_rects.append(rect)
                if len(containing_rects) == 2:
                    overlaps.add(tuple(sorted(containing_rects)))

    total_rectangles = len(rectangles) + len(overlaps)
    print(total_rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

find_rectangles(grid)