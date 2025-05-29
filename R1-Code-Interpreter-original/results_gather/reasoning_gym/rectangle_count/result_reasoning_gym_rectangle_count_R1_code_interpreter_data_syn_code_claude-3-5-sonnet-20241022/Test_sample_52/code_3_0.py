def count_rectangles(grid):
    # Convert grid to list of lines
    lines = grid.strip('\n').split('\n')
    
    # Find all potential corners (# or █)
    corners = []
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            if lines[i][j] in '#█':
                corners.append((i, j))
    
    def is_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if (top, left) not in corners or (top, right) not in corners or \
           (bottom, left) not in corners or (bottom, right) not in corners:
            return False
        
        # Check horizontal edges
        for j in range(left + 1, right):
            if lines[top][j] not in '#█' or lines[bottom][j] not in '#█':
                return False
        
        # Check vertical edges
        for i in range(top + 1, bottom):
            if lines[i][left] not in '#█' or lines[i][right] not in '#█':
                return False
        
        # Check if it's a minimal rectangle (no other corners inside)
        for i, j in corners:
            if top < i < bottom and left < j < right:
                return False
        
        return True
    
    # Find all valid rectangles
    rectangles = set()
    for i1, j1 in corners:
        for i2, j2 in corners:
            if i2 > i1 and j2 > j1:  # potential top-left and bottom-right corners
                if is_rectangle(i1, j1, i2, j2):
                    rectangles.add((i1, j1, i2, j2))
    
    print(len(rectangles))

# Your ASCII grid
grid = """    #################################################################           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                              #################█########   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                   ###########████#############█########   
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             ########    
    ####################################█#############█###############     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        ###############              #     #    
                                                                     #######    
                                                                      ########  
                                                                      #      #  
                                                                      ########  """

count_rectangles(grid)