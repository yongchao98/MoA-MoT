def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Function to check if a point forms top-left corner of a rectangle
    def is_rectangle(r, c):
        if grid[r][c] != '#':
            return False
            
        # Find width (move right until we find another '#')
        width = 1
        while c + width < cols and grid[r][c + width] == '#':
            width += 1
        
        # Find height (move down until we find another '#')
        height = 1
        while r + height < rows and grid[r + height][c] == '#':
            height += 1
            
        # Verify if it forms a complete rectangle
        if width < 2 or height < 2:
            return False
            
        # Check all corners and edges
        for i in range(r, r + height):
            if grid[i][c] != '#' or grid[i][c + width - 1] != '#':
                return False
                
        for j in range(c, c + width):
            if grid[r][j] != '#' or grid[r + height - 1][j] != '#':
                return False
                
        return True
    
    # Scan the grid for rectangles
    for i in range(rows):
        for j in range(cols):
            if is_rectangle(i, j):
                count += 1
    
    return count

# Create the grid as a list of strings
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                   ####                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   ####                                         ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))