def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_corner(i, j):
        if i < 0 or j < 0 or i >= height or j >= width:
            return False
        return grid[i][j] in '#█'
    
    def check_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if not all(is_valid_corner(y, x) for y, x in [(top, left), (top, right), (bottom, left), (bottom, right)]):
            return False
        
        # Check horizontal borders
        for x in range(left + 1, right):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        
        # Check vertical borders
        for y in range(top + 1, bottom):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        
        # Check if there's any content inside (must be empty or have another rectangle)
        has_interior = False
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] not in ' #█':
                    return False
                if grid[y][x] in '#█':
                    has_interior = True
        
        # For very small rectangles (2x2), we don't need interior
        if (bottom - top > 2 and right - left > 2) and not has_interior:
            return False
            
        return True
    
    def is_minimal_rectangle(top, left, bottom, right):
        # Check if this rectangle can be split into smaller valid rectangles
        for y in range(top + 1, bottom):
            all_border = True
            for x in range(left + 1, right):
                if grid[y][x] not in '#█':
                    all_border = False
                    break
            if all_border:
                return False
                
        for x in range(left + 1, right):
            all_border = True
            for y in range(top + 1, bottom):
                if grid[y][x] not in '#█':
                    all_border = False
                    break
            if all_border:
                return False
        return True
    
    rectangles = set()
    # Find all potential rectangles
    for top in range(height):
        for left in range(width):
            if not is_valid_corner(top, left):
                continue
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if check_rectangle(top, left, bottom, right) and is_minimal_rectangle(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Input grid (same as before)
grid = [
    # ... (same grid content as before)
]

print(find_rectangles(grid))