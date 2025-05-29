def is_corner(grid, i, j):
    # Check if a position is a corner (has perpendicular edges)
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    
    char = grid[i][j]
    if char not in ['#', '█']:
        return False
        
    # Count adjacent edges
    edges = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for di, dj in directions:
        ni, nj = i + di, j + dj
        if (0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and 
            grid[ni][nj] in ['#', '█']):
            edges += 1
    
    return edges >= 2

def find_rectangle(grid, start_i, start_j):
    # Find complete rectangle starting from top-left corner
    if not is_corner(grid, start_i, start_j):
        return None
    
    # Find right edge
    right_j = start_j
    while right_j < len(grid[0]) and grid[start_i][right_j] in ['#', '█']:
        right_j += 1
    right_j -= 1
    
    # Find bottom edge
    bottom_i = start_i
    while bottom_i < len(grid) and grid[bottom_i][start_j] in ['#', '█']:
        bottom_i += 1
    bottom_i -= 1
    
    # Verify bottom-right corner
    if not is_corner(grid, bottom_i, right_j):
        return None
    
    # Verify complete rectangle
    for i in range(start_i, bottom_i + 1):
        if grid[i][start_j] not in ['#', '█'] or grid[i][right_j] not in ['#', '█']:
            return None
    for j in range(start_j, right_j + 1):
        if grid[start_i][j] not in ['#', '█'] or grid[bottom_i][j] not in ['#', '█']:
            return None
    
    return (start_i, start_j, bottom_i, right_j)

def count_rectangles(grid):
    rectangles = set()
    grid = [list(row) for row in grid.split('\n')]
    
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['#', '█']:
                rect = find_rectangle(grid, i, j)
                if rect:
                    rectangles.add(rect)
    
    return len(rectangles)

# The ASCII grid
grid = """    #################################################                           
    #                                               #                           
    #################################################                           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                        ####### 
        ################################################################█#    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
        #                                                               ##    # 
       #█#############################################################  ##    # 
       ##                                                            #  ##    # 
       ##                                        ####################█##██##  # 
       ##                                        #                   #  ## #  # 
       ##                                        #                   #  ## #  # 
       ##                                        #                   #  ## #  # 
       ##                                        #                   #  ## #  # 
       ##                                        #                   #  #█#█### 
       ##                                        #                   #   # #    
       ##                                        #                   #   █#█    
       ##                                        #                   #   █ █    
       ##                                        #                   #   █ █    
       ##                                        #                   #   █ █    
       ##                                        #                   #   █ █    
       ##                                        #                   #   █ █    
       ##                                        #                   #   █ █    
       ##########################################█###################█###█ █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         #                   #   # █    
       #                                         ####################█###█#█    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       #                                                             #   # #    
       ###############################################################   ###    
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(count_rectangles(grid))