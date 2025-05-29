def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_char(i, j):
        if i < 0 or i >= height or j < 0 or j >= width:
            return False
        return grid[i][j] in '#█'
    
    def find_rectangle_corners(i, j):
        # Find right edge
        right = j
        while right < width and is_valid_char(i, right):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = i
        while bottom < height and is_valid_char(bottom, j):
            bottom += 1
        bottom -= 1
        
        # Verify it's a complete rectangle
        for x in range(i, bottom + 1):
            if not is_valid_char(x, right):
                return None
        for y in range(j, right + 1):
            if not is_valid_char(bottom, y):
                return None
            
        return (i, j, bottom, right)
    
    rectangles = set()
    overlap_count = 0
    
    # Find all rectangles
    for i in range(height):
        for j in range(width):
            if grid[i][j] in '#█':
                # Check if this could be top-left corner
                if (not is_valid_char(i-1, j)) and (not is_valid_char(i, j-1)):
                    rect = find_rectangle_corners(i, j)
                    if rect:
                        rectangles.add(rect)
    
    # Count overlapping points
    for i in range(height):
        for j in range(width):
            if grid[i][j] == '█':
                overlap_count += 1
    
    # Each overlapping point (█) represents where two rectangles intersect
    # So we need to add one extra rectangle for each overlap
    total = len(rectangles) + overlap_count
    print(total)

# The ASCII grid (same as before)
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