def is_valid_rectangle(grid, x1, y1, x2, y2):
    # Check if the corners are marked
    corners = [(x1, y1), (x1, y2), (x2, y1), (x2, y2)]
    for x, y in corners:
        if grid[y][x] not in '#█':
            return False
    
    # Check horizontal edges
    for y in [y1, y2]:
        for x in range(x1 + 1, x2):
            if grid[y][x] not in '#█':
                return False
    
    # Check vertical edges
    for x in [x1, x2]:
        for y in range(y1 + 1, y2):
            if grid[y][x] not in '#█':
                return False
    
    return True

def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert input to list of strings for easier processing
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] not in '#█':
                continue
            
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if grid[y2][x2] not in '#█':
                        continue
                        
                    if is_valid_rectangle(grid, x1, y1, x2, y2):
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

print(count_rectangles(grid))