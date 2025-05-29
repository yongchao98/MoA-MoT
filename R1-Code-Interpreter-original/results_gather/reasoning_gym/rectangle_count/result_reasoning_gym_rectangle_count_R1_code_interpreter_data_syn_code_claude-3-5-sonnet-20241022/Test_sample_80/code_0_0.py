def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    visited = set()

    def is_corner(i, j):
        if i < 0 or i >= height or j < 0 or j >= width:
            return False
        char = grid[i][j]
        return char in '#█'

    def trace_rectangle(start_i, start_j):
        # Find right edge
        right_j = start_j
        while right_j < width and grid[start_i][right_j] in '#█':
            right_j += 1
        right_j -= 1

        # Find bottom edge
        bottom_i = start_i
        while bottom_i < height and grid[bottom_i][start_j] in '#█':
            bottom_i += 1
        bottom_i -= 1

        # Verify rectangle
        for i in range(start_i, bottom_i + 1):
            if grid[i][start_j] not in '#█' or grid[i][right_j] not in '#█':
                return None
        for j in range(start_j, right_j + 1):
            if grid[start_i][j] not in '#█' or grid[bottom_i][j] not in '#█':
                return None

        return (start_i, start_j, bottom_i, right_j)

    def add_rectangle(rect):
        if rect is None:
            return
        start_i, start_j, bottom_i, right_j = rect
        
        # Check if this rectangle is unique
        for existing in rectangles:
            if rect == existing:
                return
            
        rectangles.append(rect)

    # Scan for potential top-left corners
    for i in range(height):
        for j in range(width):
            if grid[i][j] in '#█':
                # Check if it's a top-left corner
                if (i == 0 or grid[i-1][j] not in '#█') and (j == 0 or grid[i][j-1] not in '#█'):
                    rect = trace_rectangle(i, j)
                    add_rectangle(rect)

                # For overlapping areas (█), check additional rectangles
                if grid[i][j] == '█':
                    # Check additional possible rectangles
                    if i > 0 and j > 0 and grid[i-1][j-1] in '#█':
                        rect = trace_rectangle(i-1, j-1)
                        add_rectangle(rect)
                    if i > 0:
                        rect = trace_rectangle(i-1, j)
                        add_rectangle(rect)
                    if j > 0:
                        rect = trace_rectangle(i, j-1)
                        add_rectangle(rect)

    return len(rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")