def is_valid_rectangle(grid, x1, y1, x2, y2):
    # Check if corners are marked
    corners = [(x1, y1), (x1, y2), (x2, y1), (x2, y2)]
    for x, y in corners:
        if grid[y][x] not in ['#', '█']:
            return False
    
    # Check horizontal edges
    for x in range(x1 + 1, x2):
        if grid[y1][x] not in ['#', '█'] and grid[y2][x] not in ['#', '█']:
            return False
    
    # Check vertical edges
    for y in range(y1 + 1, y2):
        if grid[y][x1] not in ['#', '█'] and grid[y][x2] not in ['#', '█']:
            return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all possible corners
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in ['#', '█']:
                corners.append((x, y))
    
    # Check all possible rectangle combinations
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            if x2 > x1 and y2 > y1:  # Ensure correct orientation
                if is_valid_rectangle(grid, x1, y1, x2, y2):
                    rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

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