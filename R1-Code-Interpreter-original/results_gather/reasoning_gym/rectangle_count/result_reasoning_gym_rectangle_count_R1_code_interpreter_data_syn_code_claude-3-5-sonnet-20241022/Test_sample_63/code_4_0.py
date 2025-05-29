def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Store only the coordinates of '#' and '█' characters
    borders = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in '#█':
                borders.add((i, j))
    
    # Function to check if a point is part of border
    def is_valid_point(point):
        return point in borders
    
    # Check each point as potential top-left corner
    for top in range(rows):
        for left in range(cols):
            if not is_valid_point((top, left)):
                continue
                
            # Skip if not likely to be top-left corner
            if is_valid_point((top-1, left)) or is_valid_point((top, left-1)):
                continue
            
            # Find right edge
            right = left + 1
            while right < cols and is_valid_point((top, right)):
                right += 1
            right -= 1
            
            if right <= left:
                continue
                
            # Find bottom edge
            bottom = top + 1
            while bottom < rows and is_valid_point((bottom, left)):
                bottom += 1
            bottom -= 1
            
            if bottom <= top:
                continue
            
            # Quick check of corners
            if not all(is_valid_point(p) for p in [
                (top, left), (top, right),
                (bottom, left), (bottom, right)
            ]):
                continue
            
            # Verify edges
            valid = True
            for x in range(left, right + 1):
                if not (is_valid_point((top, x)) and is_valid_point((bottom, x))):
                    valid = False
                    break
            if not valid:
                continue
                
            for y in range(top, bottom + 1):
                if not (is_valid_point((y, left)) and is_valid_point((y, right))):
                    valid = False
                    break
            if valid:
                count += 1

    return count

# Read grid
grid = [line for line in iter(input, '') if line.strip()]
print(count_rectangles(grid))