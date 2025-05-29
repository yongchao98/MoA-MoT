def find_minimal_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def trace_rectangle(i, j):
        if lines[i][j] not in '#█':
            return None
            
        # Find right edge
        right = j
        while right < width and lines[i][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = i
        while bottom < height and lines[bottom][j] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        for x in range(i, bottom + 1):
            if lines[x][right] not in '#█':
                return None
        for y in range(j, right + 1):
            if lines[bottom][y] not in '#█':
                return None
                
        return (i, j, bottom, right)
    
    # Find all minimal rectangles
    rectangles = set()
    visited = set()
    
    for i in range(height):
        for j in range(width):
            if (i, j) not in visited and lines[i][j] in '#█':
                rect = trace_rectangle(i, j)
                if rect:
                    rectangles.add(rect)
                    # Mark all points in the rectangle as visited
                    for x in range(rect[0], rect[2] + 1):
                        for y in range(rect[1], rect[3] + 1):
                            visited.add((x, y))
    
    return len(rectangles)

# Test grid (same as before)
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

print(find_minimal_rectangles(grid))