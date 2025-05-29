def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Helper function to check if a point is a corner
    def is_corner(y, x):
        if y < 0 or y >= height or x < 0 or x >= width:
            return False
        return grid[y][x] in ['#', '█']
    
    # Helper function to validate a rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check corners
        if not (is_corner(top, left) and is_corner(top, right) and 
                is_corner(bottom, left) and is_corner(bottom, right)):
            return False
        
        # Check top and bottom edges
        for x in range(left + 1, right):
            if grid[top][x] not in ['#', '█']:
                return False
            if grid[bottom][x] not in ['#', '█']:
                return False
        
        # Check left and right edges
        for y in range(top + 1, bottom):
            if grid[y][left] not in ['#', '█']:
                return False
            if grid[y][right] not in ['#', '█']:
                return False
        
        return True
    
    # Find top-left corners and trace rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in ['#', '█']:
                # Look for right edge
                right = x + 1
                while right < width and grid[y][right] in ['#', '█']:
                    right += 1
                right -= 1
                
                # Look for bottom edge
                bottom = y + 1
                while bottom < height and grid[bottom][x] in ['#', '█']:
                    bottom += 1
                bottom -= 1
                
                # Validate rectangle
                if validate_rectangle(y, x, bottom, right):
                    # Check if this is likely a new rectangle (not counted before)
                    # by checking if it's a top-left corner
                    if (y == 0 or grid[y-1][x] not in ['#', '█']) and \
                       (x == 0 or grid[y][x-1] not in ['#', '█']):
                        count += 1
    
    return count

# Parse the input grid
grid_str = """                                            #########                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #                           
                                            #       #             ##############
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                       #####█#######█#############█########    #
                                       #    #       #             #       #    #
                                       #    #       #             #       #    #
                                       #    #       #             #       #    #
                                       #    #       #             #       #    #
                                       #####█#######█#############█########    #
                                            #       #             #            #
                                            #       #             #            #
                                            #       #             #            #
                                            #########             #            #
                                                                  ##############
                                                                  ##          ##
                    ##############################################██######### ##
                    #                    ##########               ##        # ##
                    #                    #        #               ##        # ##
                    #                    #        #               ##        # ##
                    #                    #        #               ##        # ##
                    #                    #        #               ##        # ##
                  ##█######              #        #               ##        # ##
                  # #     #              #        #               ##        # ##
                  # #     #              #        #               ##        # ##
                  # #     #              #        #               #█########█#█#
                  # #     #              #        #                #        # # 
                  # ######█##############█########█################█########█## 
                  # #█####█##############█########█################█########### 
                  #  #    #              #        #                #         ## 
                  #  #    #              ##########                #         ## 
                  #  #    #                                        #         ## 
                  #  #    #                                        #         ## 
                  #  #    #                                        #         ## 
                  #  #    #                                        ##########█# 
                  ###█#####                                                  #  
                     #                            ####################       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            #                  #       #  
                     #                            ############████████########  
                     #########################################█#############█#  
                                                              #             #   
                                                              ###############   
                                                                                """

# Convert the grid string to a 2D list
grid = [list(line) for line in grid_str.split('\n')]

# Count rectangles
count = find_rectangles(grid)
print(count)