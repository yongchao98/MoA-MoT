def is_valid_rectangle(grid, x1, y1, x2, y2):
    # Check if points form a valid rectangle by verifying borders
    if x1 >= x2 or y1 >= y2:
        return False
    
    # Check top and bottom borders
    for x in range(x1, x2 + 1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
    
    # Check left and right borders
    for y in range(y1, y2 + 1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible corners
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corners.append((x, y))
    
    # Check all possible rectangle combinations
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i:]:
            if x2 >= x1 and y2 >= y1:
                if is_valid_rectangle(grid, x1, y1, x2, y2):
                    rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# The ASCII grid (your input)
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                    ####################################                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #                        
                    #                                  #    ####                
         ###########█##################################█####█##█#######         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    #  #      #         
         #          #                                  #    ####      #         
         #          #                                  #              #         
         #          #                                  #              #         
         #     #####█#####                             #              #         
         #     #    #    #                             #              #         
         #     #    #    #                             #              #         
         ######█####█####█#############################█###############         
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #                        
               #    #    #                             #              ########  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              #      #  
               #    #    #                             #              ########  
               #    #    #          ###################█#                       
               #    #    #          #                  ##                       
               #    #    #          #                  #█################       
               #    #    #          #                  #█               #       
               #####█#####          #                  #█               #       
                    ################█████████████████████               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #               #       
                                                        #################       
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))