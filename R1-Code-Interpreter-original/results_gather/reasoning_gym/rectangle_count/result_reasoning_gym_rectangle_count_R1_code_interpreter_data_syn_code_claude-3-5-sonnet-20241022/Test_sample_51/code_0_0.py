def is_rectangle(grid, x, y, width, height):
    # Check if the given dimensions form a valid rectangle
    for i in range(x, x + width):
        if grid[y][i] not in '#█' or grid[y + height - 1][i] not in '#█':
            return False
    for i in range(y, y + height):
        if grid[i][x] not in '#█' or grid[i][x + width - 1] not in '#█':
            return False
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Try all possible widths and heights
                for w in range(2, width - x + 1):
                    for h in range(2, height - y + 1):
                        if x + w <= width and y + h <= height:
                            if is_rectangle(grid, x, y, w, h):
                                # Store rectangle coordinates
                                rectangles.add((x, y, w, h))
    
    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                      ######################                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #              #####                 
                      #                    #              #   #                 
                      ######################              #   #                 
                                                          #   #                 
                                                          #   #                 
                             #############################█#  #                 
                             #                            ##  #                 
                             #                   #########██##█#######          
                             #                   #        ##  #      #          
                             ####################█########█#  #      # #########
                                                 #        #   #      # #       #
    #########                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #########                                    #        #   #      # #       #
                                                 #        #####      # #       #
                                                 #                   # #       #
                                                 #                   # #       #
                                                 #  #######          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #########
                                                 #  #     #          #          
   #####                                         #  #     #          #          
   #   #                                         ###█#####█###########          
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #######                     
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #####                                                                        
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))