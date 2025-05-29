def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    def is_valid_edge_char(char):
        return char in ['#', '█']
    
    def trace_rectangle(start_y, start_x):
        # Skip if not a valid starting point
        if not is_valid_edge_char(grid[start_y][start_x]):
            return False
            
        # Find right edge
        right_x = start_x
        while right_x < width and is_valid_edge_char(grid[start_y][right_x]):
            right_x += 1
        right_x -= 1
        
        # Find bottom edge
        bottom_y = start_y
        while bottom_y < height and is_valid_edge_char(grid[bottom_y][start_x]):
            bottom_y += 1
        bottom_y -= 1
        
        # Validate rectangle
        # Check if all corners exist
        if not (is_valid_edge_char(grid[start_y][right_x]) and 
                is_valid_edge_char(grid[bottom_y][start_x]) and 
                is_valid_edge_char(grid[bottom_y][right_x])):
            return False
            
        # Check all edges
        # Top edge
        for x in range(start_x, right_x + 1):
            if not is_valid_edge_char(grid[start_y][x]):
                return False
                
        # Bottom edge
        for x in range(start_x, right_x + 1):
            if not is_valid_edge_char(grid[bottom_y][x]):
                return False
                
        # Left edge
        for y in range(start_y, bottom_y + 1):
            if not is_valid_edge_char(grid[y][start_x]):
                return False
                
        # Right edge
        for y in range(start_y, bottom_y + 1):
            if not is_valid_edge_char(grid[y][right_x]):
                return False
                
        # Check if this is a new rectangle by verifying it's a top-left corner
        if (start_y == 0 or not is_valid_edge_char(grid[start_y-1][start_x])) and \
           (start_x == 0 or not is_valid_edge_char(grid[start_y][start_x-1])):
            return True
            
        return False
    
    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            if trace_rectangle(y, x):
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