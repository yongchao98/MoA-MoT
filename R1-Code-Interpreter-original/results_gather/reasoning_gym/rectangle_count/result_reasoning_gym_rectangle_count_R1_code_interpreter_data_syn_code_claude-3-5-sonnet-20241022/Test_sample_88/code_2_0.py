def is_complete_rectangle(grid, x1, y1, x2, y2):
    # Validate dimensions
    if x2 <= x1 or y2 <= y1:
        return False
        
    # Check top and bottom edges
    for x in range(x1, x2 + 1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
            
    # Check left and right edges
    for y in range(y1, y2 + 1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
            
    # Check if there's actual content inside (not just a line)
    if x2 - x1 < 1 or y2 - y1 < 1:
        return False
        
    return True

def find_rectangles(grid):
    # Clean up the grid
    lines = [line.rstrip() for line in grid.split('\n') if line.strip()]
    grid = [list(line) for line in lines]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all corner points
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corners.append((x, y))
    
    # Try all possible corner combinations
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            if x2 > x1 and y2 > y1:  # Ensure proper rectangle orientation
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