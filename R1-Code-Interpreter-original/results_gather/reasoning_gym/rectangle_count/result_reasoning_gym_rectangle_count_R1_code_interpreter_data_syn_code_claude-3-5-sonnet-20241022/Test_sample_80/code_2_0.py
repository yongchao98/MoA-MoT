def count_rectangles(grid):
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    
    # Convert grid to list of strings for faster access
    grid = [list(row) for row in grid]
    
    # Find top-left corners and process rectangles from there
    count = 0
    for i in range(height):
        for j in range(width):
            # Skip if not a potential corner
            if grid[i][j] not in '#█':
                continue
            
            # Skip if not a top-left corner
            if i > 0 and grid[i-1][j] in '#█':
                continue
            if j > 0 and grid[i][j-1] in '#█':
                continue
            
            # Find width of potential rectangle
            width_rect = 0
            while j + width_rect < width and grid[i][j + width_rect] in '#█':
                width_rect += 1
            if width_rect < 2:  # Skip if too narrow
                continue
                
            # Find height of potential rectangle
            height_rect = 0
            while i + height_rect < height and grid[i + height_rect][j] in '#█':
                height_rect += 1
            if height_rect < 2:  # Skip if too short
                continue
            
            # Quick validate rectangle
            is_valid = True
            has_overlap = False
            
            # Check bottom edge
            for x in range(j, j + width_rect):
                if grid[i + height_rect - 1][x] not in '#█':
                    is_valid = False
                    break
                if grid[i + height_rect - 1][x] == '█':
                    has_overlap = True
            
            if not is_valid:
                continue
                
            # Check right edge
            for y in range(i, i + height_rect):
                if grid[y][j + width_rect - 1] not in '#█':
                    is_valid = False
                    break
                if grid[y][j + width_rect - 1] == '█':
                    has_overlap = True
            
            if not is_valid:
                continue
            
            # If we haven't found overlap yet, check interior
            if not has_overlap:
                for y in range(i, i + height_rect):
                    for x in range(j, j + width_rect):
                        if grid[y][x] == '█':
                            has_overlap = True
                            break
                    if has_overlap:
                        break
            
            # Count rectangles
            count += 1
            if has_overlap:
                count += 1
    
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        grid.append(line)
except EOFError:
    pass

# Print result
print(f"<<<{count_rectangles(grid)}>>>")