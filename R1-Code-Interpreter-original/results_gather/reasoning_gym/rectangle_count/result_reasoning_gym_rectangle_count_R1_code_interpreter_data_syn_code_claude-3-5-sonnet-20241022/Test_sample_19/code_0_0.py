def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = 0
    visited = set()

    def is_valid_char(char):
        return char in '#█'

    def trace_rectangle(start_y, start_x):
        if (start_y, start_x) in visited:
            return False

        # Find right boundary
        right_x = start_x
        while right_x < width and is_valid_char(grid[start_y][right_x]):
            right_x += 1
        right_x -= 1

        # Find bottom boundary
        bottom_y = start_y
        while bottom_y < height and is_valid_char(grid[bottom_y][start_x]):
            bottom_y += 1
        bottom_y -= 1

        # Verify rectangle
        for y in range(start_y, bottom_y + 1):
            for x in range(start_x, right_x + 1):
                if not is_valid_char(grid[y][x]):
                    return False
                if x in (start_x, right_x) or y in (start_y, bottom_y):
                    if grid[y][x] not in '#█':
                        return False

        # Mark as visited
        for y in range(start_y, bottom_y + 1):
            for x in range(start_x, right_x + 1):
                if grid[y][x] == '#':
                    visited.add((y, x))

        return True

    # Convert input string to grid
    grid = [list(line) for line in grid.split('\n')]

    # Find all potential rectangles
    for y in range(height):
        for x in range(width):
            if is_valid_char(grid[y][x]) and grid[y][x] == '#':
                if trace_rectangle(y, x):
                    rectangles += 1

    # Count overlapping rectangles (█)
    overlap_points = [(y, x) for y in range(height) for x in range(width) if grid[y][x] == '█']
    for y, x in overlap_points:
        if (y, x) not in visited:
            rectangles += 2  # Each █ represents two overlapping rectangles

    print(rectangles)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                ###########     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                                                                #         #     
                             ###################################█##       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  # #       #     
                             #                                  ##█########     
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #             
                             #                                    #        #####
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             #                                    #        #   #
                             ######################################        #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                           #   #
                                                                      #####█###█
                                                                      #    #   █
                                                                      #    ####█
                                                                      #        #
                                                                      #        #
                                                                      #        #
                                                                      #        #
                                                                      ##########
                                                                                
                                                                                """

find_rectangles(grid)