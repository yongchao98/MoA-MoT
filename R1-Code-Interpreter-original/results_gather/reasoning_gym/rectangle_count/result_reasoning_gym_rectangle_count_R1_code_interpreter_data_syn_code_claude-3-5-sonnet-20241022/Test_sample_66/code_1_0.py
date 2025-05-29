def is_valid_char(c):
    return c in ['#', '█']

def find_rectangles(grid):
    lines = grid.split('\n')
    height = len(lines)
    width = len(lines[0])
    count = 0
    
    # Convert to 2D array for easier processing
    grid_array = [list(line) for line in lines]
    
    # Function to check if a rectangle exists at given coordinates
    def is_complete_rectangle(top, left, bottom, right):
        # Check if all corners exist
        corners = [(top, left), (top, right), (bottom, left), (bottom, right)]
        for i, j in corners:
            if not is_valid_char(grid_array[i][j]):
                return False
        
        # Check horizontal edges
        for j in range(left, right + 1):
            if not is_valid_char(grid_array[top][j]) or not is_valid_char(grid_array[bottom][j]):
                return False
                
        # Check vertical edges
        for i in range(top, bottom + 1):
            if not is_valid_char(grid_array[i][left]) or not is_valid_char(grid_array[i][right]):
                return False
                
        # Check if interior is mostly empty (except for overlaps)
        has_interior = False
        for i in range(top + 1, bottom):
            for j in range(left + 1, right):
                if grid_array[i][j] not in [' ', '█']:
                    return False
                if grid_array[i][j] == ' ':
                    has_interior = True
        
        return has_interior

    # Find all possible rectangles
    rectangles = set()
    for top in range(height):
        for left in range(width):
            if not is_valid_char(grid_array[top][left]):
                continue
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if is_complete_rectangle(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# The ASCII grid (same as before)
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

print(find_rectangles(grid))