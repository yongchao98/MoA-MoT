def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check if all corners are present
        if not all(grid[y][x] in ['#', '█'] for y, x in [
            (top, left), (top, right), (bottom, left), (bottom, right)
        ]):
            return False
            
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
        
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
            
        # Check if it's minimal (no smaller rectangle can form this)
        for x in range(left + 1, right):
            if all(grid[y][x] not in ['#', '█'] for y in range(top + 1, bottom)):
                return False
        for y in range(top + 1, bottom):
            if all(grid[y][x] not in ['#', '█'] for x in range(left + 1, right)):
                return False
                
        return True

    # Find all potential top-left corners
    for y in range(height-1):
        for x in range(width-1):
            if grid[y][x] not in ['#', '█']:
                continue
                
            # Only consider points that could be top-left corners
            if y > 0 and grid[y-1][x] in ['#', '█']:
                continue
            if x > 0 and grid[y][x-1] in ['#', '█']:
                continue
                
            # Find potential bottom-right corners
            for bottom in range(y + 1, height):
                for right in range(x + 1, width):
                    if grid[bottom][right] in ['#', '█']:
                        if is_valid_rectangle(y, x, bottom, right):
                            rectangles.add((y, x, bottom, right))

    return len(rectangles)

# Input grid (same as before)
grid = [
    # ... (same grid as before, all rows)
]

print(find_rectangles(grid))