def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    # Find all corners (potential rectangle corners)
    corners = []
    for i in range(height):
        for j in range(width):
            if lines[i][j] in '#█':
                corners.append((i, j))
    
    def is_rectangle(top, left, bottom, right):
        # Check if the points form a valid rectangle
        if top >= bottom or left >= right:
            return False
            
        # Check all borders
        for i in range(top, bottom + 1):
            if not (lines[i][left] in '#█' and lines[i][right] in '#█'):
                return False
        for j in range(left, right + 1):
            if not (lines[top][j] in '#█' and lines[bottom][j] in '#█'):
                return False
        return True
    
    # Find all valid rectangles
    rectangles = set()
    for i1, j1 in corners:
        for i2, j2 in corners:
            if i2 > i1 and j2 > j1:  # potential bottom-right corner
                if is_rectangle(i1, j1, i2, j2):
                    rectangles.add((i1, j1, i2, j2))
    
    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                  ############                  
                                                  #          #                  
                                                  #          #                  
                                                  #          #                  
                            ######################█##########█###               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
                            #                     #          #  #               
            ################█#######              #          #  #               
            #               #######█##############█##########█###               
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              #          #                  
            #                      #              ############                  
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                  ##########
            #                      #                                  #        #
            #                      #                                  #        #
            #                      #           #######################█###     #
            #                      #           #                      #  #     #
            #                      #           #                      #  #     #
            #                      #           #######################█###     #
            #                      #                                  #        #
            #                      #                                  ##########
            #                      #                                            
            #  ####################█####                                        
            #  #                   #   #                                        
            ###█####################   #                                        
               #                       #                                        
               #                       #                                        
               #                       #                                        
               #########################                                        
                                                                                """

print(find_rectangles(grid))