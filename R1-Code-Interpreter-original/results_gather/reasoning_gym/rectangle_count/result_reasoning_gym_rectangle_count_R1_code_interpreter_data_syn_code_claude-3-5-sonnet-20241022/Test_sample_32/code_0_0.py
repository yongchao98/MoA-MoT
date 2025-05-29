def is_border_char(c):
    return c in '#█'

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if not is_border_char(grid[y][x]):
                continue
                
            # For each potential corner, try to find complete rectangle
            # First find right edge
            right = x + 1
            while right < width and is_border_char(grid[y][right]):
                right += 1
            if right == x + 1:  # No horizontal line found
                continue
                
            # Then find bottom edge
            bottom = y + 1
            valid = True
            while bottom < height and is_border_char(grid[bottom][x]):
                # Verify right edge exists at each level
                if not is_border_char(grid[bottom][right-1]):
                    valid = False
                    break
                # Verify top and bottom edges are complete
                for ix in range(x+1, right-1):
                    if not is_border_char(grid[y][ix]) or not is_border_char(grid[bottom][ix]):
                        valid = False
                        break
                if not valid:
                    break
                bottom += 1
                
            if valid and bottom > y + 1:
                # Found a valid rectangle
                rectangles.append((x, y, right-1, bottom-1))
    
    # Process overlapping rectangles
    final_count = len(rectangles)
    
    # Check for overlapping rectangles (marked by █)
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                # Count how many rectangles this pixel belongs to
                overlap_count = sum(1 for (x1, y1, x2, y2) in rectangles 
                                  if x1 <= x <= x2 and y1 <= y <= y2)
                if overlap_count > 1:
                    # Add extra rectangle for each overlap point
                    final_count += overlap_count - 1

    return final_count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")