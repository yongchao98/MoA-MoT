def count_rectangles(grid):
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Store coordinates of all '#' and '█' characters
    points = set()
    overlaps = 0
    
    # Single pass to collect points and count overlaps
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                points.add((y, x))
                if grid[y][x] == '█':
                    overlaps += 1
    
    # Function to check if a rectangle exists
    def is_rectangle(top, left, bottom, right):
        # Check all corners exist
        if not all((top, left) in points,
                   (top, right) in points,
                   (bottom, left) in points,
                   (bottom, right) in points):
            return False
        
        # Check horizontal edges
        for x in range(left + 1, right):
            if (top, x) not in points or (bottom, x) not in points:
                return False
                
        # Check vertical edges
        for y in range(top + 1, bottom):
            if (y, left) not in points or (y, right) not in points:
                return False
        
        return True
    
    # Find rectangles by checking only valid corner pairs
    corners_y = sorted(y for y, x in points)
    corners_x = sorted(x for y, x in points)
    
    # Only check unique y and x coordinates
    unique_y = list(dict.fromkeys(corners_y))
    unique_x = list(dict.fromkeys(corners_x))
    
    # Check potential rectangles
    for i, top in enumerate(unique_y[:-1]):
        for j, left in enumerate(unique_x[:-1]):
            if (top, left) not in points:
                continue
            # Look for next valid bottom and right coordinates
            for bottom in unique_y[i+1:]:
                for right in unique_x[j+1:]:
                    if is_rectangle(top, left, bottom, right):
                        count += 1
                        break
                break
    
    return count + overlaps

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(line)
except EOFError:
    pass

result = count_rectangles(grid)
print(f"<<<{result}>>>")