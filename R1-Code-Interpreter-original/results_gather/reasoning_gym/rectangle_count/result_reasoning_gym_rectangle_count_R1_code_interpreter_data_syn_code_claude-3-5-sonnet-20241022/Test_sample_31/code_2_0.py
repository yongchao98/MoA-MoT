def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Find top-left corners and validate rectangles immediately
    for y in range(height-1):
        for x in range(width-1):
            # Skip if not a potential corner
            if grid[y][x] not in ['#', '█']:
                continue
                
            # Look for matching bottom-right corners
            for y2 in range(y+1, height):
                for x2 in range(x+1, width):
                    # Quick check for bottom-right corner
                    if grid[y2][x2] not in ['#', '█']:
                        continue
                        
                    # Check if other corners exist
                    if grid[y][x2] not in ['#', '█'] or grid[y2][x] not in ['#', '█']:
                        continue
                    
                    # Validate edges
                    valid = True
                    
                    # Check top and bottom edges
                    for i in range(x+1, x2):
                        if grid[y][i] not in ['#', '█'] or grid[y2][i] not in ['#', '█']:
                            valid = False
                            break
                    
                    if not valid:
                        continue
                    
                    # Check left and right edges
                    for i in range(y+1, y2):
                        if grid[i][x] not in ['#', '█'] or grid[i][x2] not in ['#', '█']:
                            valid = False
                            break
                    
                    if valid:
                        count += 1
                        # Skip remaining x2 values as they can't form valid rectangles
                        break
    
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(list(line))
except EOFError:
    pass

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")