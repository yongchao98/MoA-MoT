def find_rectangle_corners(grid, start_i, start_j):
    # Find the right edge
    right_j = start_j
    while right_j < len(grid[0]) and grid[start_i][right_j] in '#█':
        right_j += 1
    right_j -= 1
    
    # Find the bottom edge
    bottom_i = start_i
    while bottom_i < len(grid) and grid[bottom_i][start_j] in '#█':
        bottom_i += 1
    bottom_i -= 1
    
    # Verify rectangle
    if not grid[start_i][right_j] in '#█' or not grid[bottom_i][start_j] in '#█' or not grid[bottom_i][right_j] in '#█':
        return None
    
    # Check all edges are complete
    for i in range(start_i, bottom_i + 1):
        if not grid[i][start_j] in '#█' or not grid[i][right_j] in '#█':
            return None
    for j in range(start_j, right_j + 1):
        if not grid[start_i][j] in '#█' or not grid[bottom_i][j] in '#█':
            return None
            
    return (start_i, start_j, bottom_i, right_j)

def count_rectangles(grid):
    rectangles = set()
    visited = set()
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in '#█' and (i, j) not in visited:
                # Try to find rectangle starting at this corner
                rect = find_rectangle_corners(grid, i, j)
                if rect:
                    rectangles.add(rect)
                    # Mark the outline as visited
                    top, left, bottom, right = rect
                    for x in range(top, bottom + 1):
                        visited.add((x, left))
                        visited.add((x, right))
                    for y in range(left, right + 1):
                        visited.add((top, y))
                        visited.add((bottom, y))
    
    return len(rectangles)

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():
            continue
        grid.append(line)
except EOFError:
    pass

# Count rectangles
result = count_rectangles(grid)
print(f"<<<{result}>>>")