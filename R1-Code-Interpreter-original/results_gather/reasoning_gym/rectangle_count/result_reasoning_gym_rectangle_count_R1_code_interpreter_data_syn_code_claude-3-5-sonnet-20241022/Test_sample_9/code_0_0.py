def is_corner(grid, i, j):
    # Check if position is a corner (has perpendicular edges)
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    
    char = grid[i][j]
    if char not in ['#', '█']:
        return False
        
    # Count different directions where we have edges
    directions = 0
    # Check right
    if j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█']:
        directions += 1
    # Check left
    if j-1 >= 0 and grid[i][j-1] in ['#', '█']:
        directions += 1
    # Check down
    if i+1 < len(grid) and grid[i+1][j] in ['#', '█']:
        directions += 1
    # Check up
    if i-1 >= 0 and grid[i-1][j] in ['#', '█']:
        directions += 1
    
    return directions == 2 and (
        (i+1 < len(grid) and grid[i+1][j] in ['#', '█'] and j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█']) or
        (i+1 < len(grid) and grid[i+1][j] in ['#', '█'] and j-1 >= 0 and grid[i][j-1] in ['#', '█']) or
        (i-1 >= 0 and grid[i-1][j] in ['#', '█'] and j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█']) or
        (i-1 >= 0 and grid[i-1][j] in ['#', '█'] and j-1 >= 0 and grid[i][j-1] in ['#', '█'])
    )

def find_rectangle(grid, start_i, start_j):
    # Find complete rectangle starting from top-left corner
    width = 0
    height = 0
    
    # Find width
    j = start_j
    while j < len(grid[0]) and grid[start_i][j] in ['#', '█']:
        width += 1
        j += 1
    
    # Find height
    i = start_i
    while i < len(grid) and grid[i][start_j] in ['#', '█']:
        height += 1
        i += 1
    
    # Verify rectangle
    for i in range(start_i, start_i + height):
        for j in range(start_j, start_j + width):
            if grid[i][j] not in ['#', '█']:
                return None
    
    return (start_i, start_j, height, width)

def count_rectangles(grid):
    rectangles = set()
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all corners
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if is_corner(grid, i, j):
                rect = find_rectangle(grid, i, j)
                if rect:
                    rectangles.add(rect)
    
    return len(rectangles)

# The ASCII grid
grid = """                                                                                
                                                                                
                               #############################################    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #          #################                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                               #          #               #                #    
                   #############          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #               #                #    
                   #          ##          #################                #    
                   #          ##                                           #    
                   #          ##                                           #    
                   #          ##                                           #    
                   #          ##                                           #    
                   #############                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #                                           #    
                               #############################################    
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                         #####  
                                                                         #   #  
                                                                         #   #  
                                                 ########################█###█# 
                                #################█#######################█## ## 
               #################█################███████████████████#####█#█#█# 
               #                ###################################█#####███##  
               #####################################################            
                                                                                
                                                                                """

print(count_rectangles(grid))