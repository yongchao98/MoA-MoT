def is_valid_char(c):
    return c in '#█'

def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # For each potential top-left corner
    for top in range(height):
        for left in range(width):
            if not is_valid_char(grid[top][left]):
                continue
                
            # For each potential bottom-right corner
            for bottom in range(top, height):
                for right in range(left, width):
                    # Check if we can form a rectangle
                    valid = True
                    
                    # Check top and bottom edges
                    for x in range(left, right + 1):
                        if not is_valid_char(grid[top][x]) or not is_valid_char(grid[bottom][x]):
                            valid = False
                            break
                    
                    if not valid:
                        continue
                    
                    # Check left and right edges
                    for y in range(top, bottom + 1):
                        if not is_valid_char(grid[y][left]) or not is_valid_char(grid[y][right]):
                            valid = False
                            break
                    
                    if valid:
                        # Add rectangle coordinates to set
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

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
print(count)