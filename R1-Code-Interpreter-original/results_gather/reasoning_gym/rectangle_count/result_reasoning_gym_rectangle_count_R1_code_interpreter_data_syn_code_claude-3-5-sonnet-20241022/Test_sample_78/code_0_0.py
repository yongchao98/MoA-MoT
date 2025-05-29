def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(i, j):
        if i < 0 or i >= height or j < 0 or j >= width:
            return False
        char = grid[i][j]
        return char in '#█'
    
    def trace_rectangle(start_i, start_j):
        # Find the right edge
        right_j = start_j
        while right_j < width and is_corner(start_i, right_j):
            right_j += 1
        right_j -= 1
        
        # Find the bottom edge
        bottom_i = start_i
        while bottom_i < height and is_corner(bottom_i, start_j):
            bottom_i += 1
        bottom_i -= 1
        
        # Verify rectangle
        for i in range(start_i, bottom_i + 1):
            if not is_corner(i, right_j):
                return None
        for j in range(start_j, right_j + 1):
            if not is_corner(bottom_i, j):
                return None
        
        return (start_i, start_j, bottom_i, right_j)
    
    # Find all rectangles
    rectangles = set()
    overlaps = 0
    
    for i in range(height):
        for j in range(width):
            if grid[i][j] in '#█':
                # Check if this could be top-left corner
                if (not is_corner(i-1, j)) and (not is_corner(i, j-1)):
                    rect = trace_rectangle(i, j)
                    if rect:
                        rectangles.add(rect)
                # Count overlaps
                if grid[i][j] == '█':
                    overlaps += 1

    # Each '█' character represents an overlap of 2 rectangles
    # So we need to add one extra rectangle for each overlap
    total_rectangles = len(rectangles) + overlaps

    print(total_rectangles)

# The ASCII grid
grid = """
                   ###############################                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #                              
                   #                             #    #####################     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                             #    #                   #     
                   #                #############█####█#######            #     
                   #                #            #    #      #            #     
                   #                #            #    #######█#############     
                   #                #            #           #                  
                   #                #            #           #                  
                   #                #            #           #                  
                   #                #            #           #              ### 
                   #                #            #           #              # # 
                   #                #############█############              # # 
                   #                             #                          # # 
                   #                             #                          # # 
                   #                             #                          # # 
                   #                             #  ####################    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
                   #                             #  #                  #    # # 
     ##############█##                           #  #                  #    # # 
     #             # #                           #  #                  #    # # 
     #             # #                           #  #                  #    # # 
     ##############█##                           #  ####################    # # 
                   #               ##############█##########################██# 
                   ################█##############                           #  
                                   #                                         #  
                                   ###########################################  
"""

find_rectangles(grid)