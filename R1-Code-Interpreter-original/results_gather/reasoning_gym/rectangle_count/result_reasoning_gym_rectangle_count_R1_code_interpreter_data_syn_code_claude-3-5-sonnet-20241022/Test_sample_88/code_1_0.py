def is_valid_corner(grid, x, y):
    return 0 <= y < len(grid) and 0 <= x < len(grid[0]) and grid[y][x] in '#█'

def is_complete_rectangle(grid, x1, y1, x2, y2):
    # Ensure proper rectangle dimensions
    if x2 <= x1 + 1 or y2 <= y1 + 1:
        return False
    
    # Check corners
    corners = [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]
    if not all(is_valid_corner(grid, x, y) for x, y in corners):
        return False
    
    # Check horizontal edges (top and bottom)
    for x in range(x1 + 1, x2):
        if not (grid[y1][x] in '#█' and grid[y2][x] in '#█'):
            return False
    
    # Check vertical edges (left and right)
    for y in range(y1 + 1, y2):
        if not (grid[y][x1] in '#█' and grid[y][x2] in '#█'):
            return False
    
    # Verify it's not just a line
    if x2 - x1 < 2 or y2 - y1 < 2:
        return False
    
    # Check if there's empty space inside
    has_empty = False
    for y in range(y1 + 1, y2):
        for x in range(x1 + 1, x2):
            if grid[y][x] == ' ':
                has_empty = True
                break
        if has_empty:
            break
    
    return has_empty

def find_rectangles(grid):
    # Clean and prepare grid
    lines = [line for line in grid.split('\n') if line.strip()]
    grid = [list(line) for line in lines]
    height = len(grid)
    width = len(grid[0])
    
    rectangles = set()
    
    # Find potential corner points
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corners.append((x, y))
    
    # Check all possible corner combinations
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            # Ensure proper orientation
            if x1 < x2 and y1 < y2:
                if is_complete_rectangle(grid, x1, y1, x2, y2):
                    rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Input grid
grid = """   ###############                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #          ###█#####################################################         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          ###█##█████████##############################█████#######         
   #             #  #       █##################            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #####                
   #             #  #       █                 #                                 
   #             #  #       █                 #                                 
   #             #  #       █                 #                                 
   ###############  #       █                 #                                 
                    #       █                 #                                 
                    ########█                 #                                 
                            #                 #                                 
                            #                 #                                 
                           #█#############    #                                 
                           ##            #    #                                 
                           ##############█#####                                 
                           #             #                                      
                           #             #                                      
                           #             #       #############################  
                           #             #       #                           #  
                           #             #       #                           #  
                           #             #       #                           #  
                           #             #       #                           #  
                           ###############       #                           #  
                               ###############   #                           #  
                               #             #   #                           #  
                               #             #   ███████████████##############  
                               #             #   #             #                
                               ###############   ###############                """

print(find_rectangles(grid))