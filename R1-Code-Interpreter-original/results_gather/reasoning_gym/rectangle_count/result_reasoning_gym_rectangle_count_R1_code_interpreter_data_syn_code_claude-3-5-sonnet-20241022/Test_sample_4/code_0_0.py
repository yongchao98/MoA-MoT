def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Helper function to check if a point is part of rectangle border
    def is_border(x, y):
        return grid[y][x] in ['#', '█']
    
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if not is_border(x, y):
                continue
                
            # Skip if this point is not a corner (needs empty space above and to the left)
            if y > 0 and is_border(x, y-1):
                continue
            if x > 0 and is_border(x-1, y):
                continue
                
            # Find width of top edge
            w = 0
            while x + w < width and is_border(x + w, y):
                w += 1
            if w < 2:  # Too small to be rectangle
                continue
                
            # Find height of left edge
            h = 0
            while y + h < height and is_border(x, y + h):
                h += 1
            if h < 2:  # Too small to be rectangle
                continue
                
            # Verify bottom and right edges exist
            valid = True
            # Check bottom edge
            for dx in range(w):
                if not is_border(x + dx, y + h - 1):
                    valid = False
                    break
            # Check right edge
            if valid:
                for dy in range(h):
                    if not is_border(x + w - 1, y + dy):
                        valid = False
                        break
                        
            if valid:
                rectangles.append((x, y, w, h))
    
    # Count overlapping rectangles
    overlaps = set()
    for i, r1 in enumerate(rectangles):
        for j, r2 in enumerate(rectangles[i+1:], i+1):
            x1, y1, w1, h1 = r1
            x2, y2, w2, h2 = r2
            
            # Check if rectangles overlap
            if (x1 < x2 + w2 and x1 + w1 > x2 and
                y1 < y2 + h2 and y1 + h1 > y2):
                # Verify overlap contains '█'
                for y in range(max(y1, y2), min(y1 + h1, y2 + h2)):
                    for x in range(max(x1, x2), min(x1 + w1, x2 + w2)):
                        if grid[y][x] == '█':
                            overlaps.add((i, j))
                            break
                    if (i, j) in overlaps:
                        break
    
    # Each overlap represents two rectangles that were counted separately
    return len(rectangles)

# Process the input grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

print(find_rectangles(grid))