def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_complete_rectangle(top, left, bottom, right):
        # Check if we have a valid rectangle border
        # Check horizontal borders
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        # Check vertical borders
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        
        # Check that inside is empty (except for possible overlapping parts)
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] == '#':
                    return False
        
        return True

    def is_minimal_rectangle(top, left, bottom, right):
        # Check if this rectangle can't be split into smaller valid rectangles
        for y in range(top + 1, bottom):
            complete_horizontal = True
            for x in range(left, right + 1):
                if grid[y][x] not in '#█':
                    complete_horizontal = False
                    break
            if complete_horizontal:
                return False
                
        for x in range(left + 1, right):
            complete_vertical = True
            for y in range(top, bottom + 1):
                if grid[y][x] not in '#█':
                    complete_vertical = False
                    break
            if complete_vertical:
                return False
        
        return True

    rectangles = set()
    # Find top-left corners
    for top in range(height):
        for left in range(width):
            if grid[top][left] in '#█':
                # Find possible bottom-right corners
                for bottom in range(top + 1, height):
                    for right in range(left + 1, width):
                        if grid[bottom][right] in '#█':
                            if is_complete_rectangle(top, left, bottom, right):
                                if is_minimal_rectangle(top, left, bottom, right):
                                    rectangles.add((top, left, bottom, right))

    return len(rectangles)

# Create the grid (same as before)
grid = [
    # ... (your full grid here with all the lines)
]

print(find_rectangles(grid))