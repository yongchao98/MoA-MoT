def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Find all corners
    corners = []
    for i in range(height-1):
        for j in range(width-1):
            if grid[i][j] in '#█':
                # Check if it's a corner (no marker above and to the left)
                if (i == 0 or grid[i-1][j] not in '#█') and (j == 0 or grid[i][j-1] not in '#█'):
                    corners.append((i, j))
    
    # Process each corner
    for corner_i, corner_j in corners:
        # Find width
        width_rect = 1
        while corner_j + width_rect < width and grid[corner_i][corner_j + width_rect] in '#█':
            width_rect += 1
        
        if width_rect < 2:
            continue
            
        # Find height
        height_rect = 1
        while corner_i + height_rect < height and grid[corner_i + height_rect][corner_j] in '#█':
            height_rect += 1
            
        if height_rect < 2:
            continue
        
        # Quick validate bottom and right edges
        valid = True
        has_overlap = False
        
        # Check bottom edge
        bottom_i = corner_i + height_rect - 1
        for j in range(corner_j, corner_j + width_rect):
            if grid[bottom_i][j] not in '#█':
                valid = False
                break
            if grid[bottom_i][j] == '█':
                has_overlap = True
                
        if not valid:
            continue
            
        # Check right edge
        right_j = corner_j + width_rect - 1
        for i in range(corner_i, corner_i + height_rect):
            if grid[i][right_j] not in '#█':
                valid = False
                break
            if grid[i][right_j] == '█':
                has_overlap = True
                
        if not valid:
            continue
            
        # If no overlap found yet, quick check for interior overlap
        if not has_overlap:
            for i in range(corner_i, corner_i + height_rect):
                if '█' in grid[i][corner_j:corner_j + width_rect]:
                    has_overlap = True
                    break
        
        # Count rectangle(s)
        count += 1 + (1 if has_overlap else 0)
    
    return count

# Read grid
grid = []
for _ in range(80):  # Read fixed number of lines
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Print result
print(f"<<<{count_rectangles(grid)}>>>")