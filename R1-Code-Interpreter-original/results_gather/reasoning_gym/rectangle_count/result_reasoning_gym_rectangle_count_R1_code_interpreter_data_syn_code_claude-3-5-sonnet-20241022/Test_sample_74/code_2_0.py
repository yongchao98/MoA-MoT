def count_rectangles(grid):
    # Convert input to list for faster processing
    grid = [list(row) for row in grid]
    rows, cols = len(grid), len(grid[0])
    count = 0
    
    # Find top-left corners and process rectangles directly
    for r in range(rows-1):
        for c in range(cols-1):
            # Skip if not a potential corner
            if grid[r][c] not in '#█':
                continue
                
            # Find width of top edge
            width = 0
            for dc in range(c, cols):
                if grid[r][dc] not in '#█':
                    break
                width += 1
            if width < 2:
                continue
                
            # Find height of left edge
            height = 0
            for dr in range(r, rows):
                if grid[dr][c] not in '#█':
                    break
                height += 1
            if height < 2:
                continue
                
            # Quick check for bottom-right corner
            if (grid[r+height-1][c+width-1] not in '#█' or
                grid[r][c+width-1] not in '#█' or
                grid[r+height-1][c] not in '#█'):
                continue
                
            # Verify rectangle edges
            valid = True
            overlap = False
            
            # Check top and bottom edges
            for dc in range(c, c+width):
                if grid[r][dc] not in '#█' or grid[r+height-1][dc] not in '#█':
                    valid = False
                    break
                if grid[r][dc] == '█' or grid[r+height-1][dc] == '█':
                    overlap = True
                    
            if not valid:
                continue
                
            # Check left and right edges
            for dr in range(r, r+height):
                if grid[dr][c] not in '#█' or grid[dr][c+width-1] not in '#█':
                    valid = False
                    break
                if grid[dr][c] == '█' or grid[dr][c+width-1] == '█':
                    overlap = True
                    
            if valid:
                count += 1
                if overlap:
                    count += 1
                    
            # Skip processed columns
            c += width-1
    
    print(count)

# Read input
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

count_rectangles(grid)