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
    # Split the grid into lines and remove empty lines
    lines = [line for line in grid.split('\n') if line]
    height = len(lines)
    width = len(lines[0]) if height > 0 else 0
    rectangles = set()
    
    # Find all possible rectangles
    for y in range(height):
        for x in range(width):
            if lines[y][x] in '#█':
                # Try all possible widths and heights
                for w in range(2, width - x + 1):
                    for h in range(2, height - y + 1):
                        if x + w <= width and y + h <= height:
                            # Check if it forms a valid rectangle
                            valid = True
                            for i in range(x, x + w):
                                if lines[y][i] not in '#█' or lines[y + h - 1][i] not in '#█':
                                    valid = False
                                    break
                            if valid:
                                for i in range(y, y + h):
                                    if lines[i][x] not in '#█' or lines[i][x + w - 1] not in '#█':
                                        valid = False
                                        break
                            if valid:
                                # Check if interior is empty (except for overlaps)
                                interior_valid = True
                                for i in range(y + 1, y + h - 1):
                                    for j in range(x + 1, x + w - 1):
                                        if lines[i][j] not in ' █':
                                            interior_valid = False
                                            break
                                    if not interior_valid:
                                        break
                                if interior_valid:
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