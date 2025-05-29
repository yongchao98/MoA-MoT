def find_rectangles(grid):
    # Convert grid to list of strings
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(y, x):
        if not (0 <= y < height and 0 <= x < width):
            return False
        if grid[y][x] not in '#█':
            return False
            
        # Check if it's a corner by looking at adjacent cells
        horizontal = 0
        vertical = 0
        
        # Check right
        if x + 1 < width and grid[y][x+1] in '#█':
            horizontal += 1
        # Check left
        if x - 1 >= 0 and grid[y][x-1] in '#█':
            horizontal += 1
        # Check down
        if y + 1 < height and grid[y+1][x] in '#█':
            vertical += 1
        # Check up
        if y - 1 >= 0 and grid[y-1][x] in '#█':
            vertical += 1
            
        return horizontal == 1 and vertical == 1

    def find_rectangle(start_y, start_x):
        # Find right edge
        right = start_x
        while right < width and grid[start_y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = start_y
        while bottom < height and grid[bottom][start_x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify all corners
        corners = [
            (start_y, start_x),
            (start_y, right),
            (bottom, start_x),
            (bottom, right)
        ]
        
        for corner_y, corner_x in corners:
            if not is_corner(corner_y, corner_x):
                return None
                
        return (start_y, start_x, bottom, right)

    rectangles = set()
    # Scan for top-left corners
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if it could be a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and \
                   (x == 0 or grid[y][x-1] not in '#█'):
                    rect = find_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)

    return len(rectangles)

# Test grid (same as before)
grid = """
                                                                                
                                                                                
                                                 ##################             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 #                #             
                                                 ##################             
                                                                                
                                                                                
                                                                                
                                                                                
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
                                                 #                     ####    #
                                                 #                     #  #    #
                                                 #                     #  #    #
                                        #########█###########          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          #  #    #
                                        #        #          #          ████#####
                                        #        #          #          #      ##
                                        #        #       ###█######    #      ##
                                        #        #       #  #     #    #      ##
                                        #        #       ###█######    #      ##
                                        #        #          #          #      ##
                                        #########█###########          #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #      ##
                                                 #                     #########
                                                 #                             #
                                                 #                             #
                                                 #        #############        #
                                                 #########█###########█#########
                                                          #           #         
                                                          #           #         
                                                          #           #         
                                                          #           #         
                   ####################███████████        #############         
                   #                   #         █                              
                   ####################███████████                              
"""

print(find_rectangles(grid))