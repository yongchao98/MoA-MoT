def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])

    def is_border(char):
        return char in '#█'

    def trace_rectangle(y, x):
        # Skip if not a potential top-left corner
        if not is_border(lines[y][x]):
            return False
        if y > 0 and is_border(lines[y-1][x]):
            return False
        if x > 0 and is_border(lines[y][x-1]):
            return False

        # Find right edge
        right = x + 1
        while right < width and is_border(lines[y][right]):
            right += 1
        right -= 1

        # Find bottom edge
        bottom = y + 1
        while bottom < height and is_border(lines[bottom][x]):
            bottom += 1
        bottom -= 1

        # Verify rectangle
        if bottom <= y or right <= x:
            return False

        # Check all borders are present
        for i in range(y, bottom + 1):
            if not is_border(lines[i][x]) or not is_border(lines[i][right]):
                return False
        for j in range(x, right + 1):
            if not is_border(lines[y][j]) or not is_border(lines[bottom][j]):
                return False

        # Check if it's a valid rectangle (has four corners)
        corners = [(y, x), (y, right), (bottom, x), (bottom, right)]
        for cy, cx in corners:
            if not is_border(lines[cy][cx]):
                return False

        return True

    count = 0
    # Scan the entire grid
    for y in range(height):
        for x in range(width):
            if trace_rectangle(y, x):
                count += 1

    return count

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