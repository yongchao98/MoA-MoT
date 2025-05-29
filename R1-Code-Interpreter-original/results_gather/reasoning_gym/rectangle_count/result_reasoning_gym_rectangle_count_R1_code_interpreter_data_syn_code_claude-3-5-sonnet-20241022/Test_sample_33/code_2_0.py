def is_complete_rectangle(grid, x1, y1, x2, y2):
    # Check if it's a minimal rectangle (no smaller rectangles inside)
    for y in range(y1 + 1, y2):
        for x in range(x1 + 1, x2):
            if grid[y][x] in '#█':
                return False
    
    # Check top and bottom borders are complete
    for x in range(x1, x2 + 1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
    
    # Check left and right borders are complete
    for y in range(y1, y2 + 1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
    
    return True

def find_minimal_rectangles(grid):
    # Split the grid into lines and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    grid = [list(line) for line in lines]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = set()
    
    # Find all border points
    borders = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                borders.append((x, y))
    
    # Check each possible rectangle
    for i, (x1, y1) in enumerate(borders):
        for j, (x2, y2) in enumerate(borders[i:], i):
            if x2 > x1 and y2 > y1:  # Ensure proper rectangle dimensions
                # Check if it forms a complete rectangle
                if is_complete_rectangle(grid, x1, y1, x2, y2):
                    # Check if it's not contained within another rectangle
                    is_minimal = True
                    for (rx1, ry1, rx2, ry2) in rectangles:
                        if (rx1 <= x1 <= rx2 and rx1 <= x2 <= rx2 and
                            ry1 <= y1 <= ry2 and ry1 <= y2 <= ry2):
                            is_minimal = False
                            break
                    if is_minimal:
                        rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Your ASCII grid input (same as before)
grid = """                                                                                
[... rest of the grid ...]
                                                                                """

print(find_minimal_rectangles(grid))