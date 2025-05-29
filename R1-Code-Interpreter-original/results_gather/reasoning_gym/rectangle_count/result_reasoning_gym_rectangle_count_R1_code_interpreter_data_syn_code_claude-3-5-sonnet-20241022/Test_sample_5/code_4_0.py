def find_rectangles(raw_grid):
    # Compress grid by removing empty columns and rows
    grid = [line for line in raw_grid if '#' in line or '█' in line]
    if not grid:
        return 0
        
    # Find all corners directly
    corners = []
    rows = len(grid)
    cols = len(grid[0])
    
    # Pre-compute valid positions
    valid_chars = {'#', '█'}
    valid_positions = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in valid_chars:
                valid_positions.add((i, j))
    
    # Find top-left corners only
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in valid_positions:
                continue
            # Check if it's a top-left corner
            is_corner = (j == 0 or grid[i][j-1] not in valid_chars) and \
                       (i == 0 or grid[i-1][j] not in valid_chars)
            if is_corner:
                corners.append((i, j))
    
    count = 0
    # For each corner, try to find the smallest valid rectangle
    for top, left in corners:
        # Find right boundary
        right = left + 1
        while right < cols and (top, right) in valid_positions:
            right += 1
        right -= 1
        
        if right == left:
            continue
            
        # Find bottom boundary
        bottom = top + 1
        while bottom < rows and (bottom, left) in valid_positions:
            bottom += 1
        bottom -= 1
        
        if bottom == top:
            continue
            
        # Quick validation of the rectangle
        valid = True
        # Check only the perimeter
        for x in range(left, right + 1):
            if (top, x) not in valid_positions or (bottom, x) not in valid_positions:
                valid = False
                break
        if valid:
            for y in range(top, bottom + 1):
                if (y, left) not in valid_positions or (y, right) not in valid_positions:
                    valid = False
                    break
        
        if valid:
            count += 1
    
    return count

# Efficient input reading
grid = []
try:
    while True:
        line = input()
        if line.strip():
            grid.append(line)
        elif grid:
            break
except EOFError:
    pass

print(f"<<<{find_rectangles(grid)}>>>")