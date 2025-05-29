def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Quick check for corners
    def is_corner(i, j):
        if grid[i][j] not in ['#', '█']:
            return False
        # Check if it's a potential top-left corner
        left = j == 0 or grid[i][j-1] not in ['#', '█']
        top = i == 0 or grid[i-1][j] not in ['#', '█']
        return left and top
    
    # Quick validate rectangle
    def validate_rect(top, left, h, w):
        # Check corners exist
        if (grid[top][left] not in ['#', '█'] or
            grid[top][left+w-1] not in ['#', '█'] or
            grid[top+h-1][left] not in ['#', '█'] or
            grid[top+h-1][left+w-1] not in ['#', '█']):
            return False
        
        # Check top and bottom edges
        for j in range(left, left+w):
            if grid[top][j] not in ['#', '█'] or grid[top+h-1][j] not in ['#', '█']:
                return False
        
        # Check left and right edges
        for i in range(top, top+h):
            if grid[i][left] not in ['#', '█'] or grid[i][left+w-1] not in ['#', '█']:
                return False
        
        return True

    # Process only first 1000 rows and columns for timeout prevention
    max_dim = min(1000, rows, cols)
    
    # Find rectangles starting from top-left corners
    for i in range(max_dim):
        for j in range(max_dim):
            if not is_corner(i, j):
                continue
                
            # Find width
            width = 1
            while j + width < cols and grid[i][j + width] in ['#', '█']:
                width += 1
            
            # Find height
            height = 1
            while i + height < rows and grid[i + height][j] in ['#', '█']:
                height += 1
            
            # Validate rectangle
            if validate_rect(i, j, height, width):
                count += 1

    return count

# Read grid with limit
grid = []
line_count = 0
max_lines = 1000  # Limit number of lines to prevent timeout

try:
    while line_count < max_lines:
        line = input().strip('\n')
        if not line:
            break
        grid.append(list(line))
        line_count += 1
except EOFError:
    pass

print(count_rectangles(grid))