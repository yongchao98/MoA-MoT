def count_rectangles(grid):
    # Convert to list and get dimensions
    grid = [list(row) for row in grid]
    H, W = len(grid), len(grid[0])
    
    # Find top-left corners and count overlaps in one pass
    corners = []  # (y, x)
    overlaps = 0
    
    for y in range(H):
        for x in range(W):
            if grid[y][x] in '#█':
                # Check if it's a top-left corner
                is_corner = (y == 0 or grid[y-1][x] not in '#█') and \
                           (x == 0 or grid[y][x-1] not in '#█')
                if is_corner:
                    corners.append((y, x))
                # Count overlaps
                if grid[y][x] == '█':
                    overlaps += 1

    # Count basic rectangles (without considering overlaps)
    rect_count = 0
    for y, x in corners:
        # Find width and height
        w = h = 0
        # Get width
        while x + w < W and grid[y][x + w] in '#█':
            w += 1
        # Get height
        while y + h < H and grid[y + h][x] in '#█':
            h += 1
        
        # Verify it's a rectangle
        if w > 1 and h > 1:
            is_rect = True
            # Quick check corners only
            if not (grid[y][x + w - 1] in '#█' and 
                   grid[y + h - 1][x] in '#█' and 
                   grid[y + h - 1][x + w - 1] in '#█'):
                is_rect = False
            
            if is_rect:
                rect_count += 1

    # Each overlap point indicates an additional rectangle
    print(rect_count + overlaps)

# Read input
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

count_rectangles(grid)