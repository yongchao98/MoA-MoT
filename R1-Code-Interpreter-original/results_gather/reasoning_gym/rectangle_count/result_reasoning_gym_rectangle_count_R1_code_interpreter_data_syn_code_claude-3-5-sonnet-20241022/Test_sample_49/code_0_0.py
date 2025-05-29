def find_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    rectangles = []

    # Function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'

    # Function to validate if coordinates are within grid
    def is_valid(y, x):
        return 0 <= y < height and 0 <= x < width

    # Function to find complete rectangle from top-left corner
    def find_rectangle(top, left):
        # Find right boundary
        right = left
        while right < width and is_border(grid[top][right]):
            right += 1
        right -= 1

        # Find bottom boundary
        bottom = top
        while bottom < height and is_border(grid[bottom][left]):
            bottom += 1
        bottom -= 1

        # Validate rectangle
        if bottom <= top or right <= left:
            return None

        # Check if all borders are present
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return None
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return None

        return (top, left, bottom, right)

    # Scan grid for rectangles
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Check if this could be top-left corner
                if (not is_valid(y-1, x) or not is_border(grid[y-1][x])) and \
                   (not is_valid(y, x-1) or not is_border(grid[y][x-1])):
                    rect = find_rectangle(y, x)
                    if rect:
                        rectangles.append(rect)

    return len(rectangles)

# Test grid
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