def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all corners (points where we have a corner of a rectangle)
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if it's a corner by looking at neighbors
                is_corner = False
                # Top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    is_corner = True
                # Top-right corner
                elif (y == 0 or grid[y-1][x] not in '#█') and (x == width-1 or grid[y][x+1] not in '#█'):
                    is_corner = True
                # Bottom-left corner
                elif (y == height-1 or grid[y+1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    is_corner = True
                # Bottom-right corner
                elif (y == height-1 or grid[y+1][x] not in '#█') and (x == width-1 or grid[y][x+1] not in '#█'):
                    is_corner = True
                
                if is_corner:
                    corners.append((y, x))
    
    # For each pair of corners, check if they form a rectangle
    for i, (y1, x1) in enumerate(corners):
        for y2, x2 in corners[i+1:]:
            if y2 > y1 and x2 > x1:  # Only check if second corner is below and to the right
                # Check if all edges exist
                valid = True
                # Check vertical edges
                for y in range(y1, y2 + 1):
                    if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                        valid = False
                        break
                if not valid:
                    continue
                
                # Check horizontal edges
                for x in range(x1, x2 + 1):
                    if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                        valid = False
                        break
                if valid:
                    rectangles.add((y1, x1, y2, x2))
    
    # Count base rectangles
    base_count = len(rectangles)
    
    # Count overlapping points
    overlap_points = sum(1 for y in range(height) for x in range(width) if grid[y][x] == '█')
    
    # Each █ represents an overlap of two rectangles, so we need to add one rectangle
    # for each overlapping point and then divide by 2 (as each overlap is counted twice)
    total_count = base_count + overlap_points // 2
    
    return total_count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
result = find_rectangles(grid)
print(result)