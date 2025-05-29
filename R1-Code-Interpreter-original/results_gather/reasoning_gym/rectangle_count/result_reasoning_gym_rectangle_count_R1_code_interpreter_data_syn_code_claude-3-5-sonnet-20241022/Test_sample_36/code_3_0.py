def find_rectangles(grid):
    if not grid:
        return 0
    
    # Convert grid to list of lists for easier handling
    grid = [list(row) for row in grid]
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    def is_valid_rect(y1, x1, y2, x2):
        # Check horizontal borders
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
        
        # Check vertical borders
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
        
        return True
    
    # Mark processed areas to avoid counting the same rectangle multiple times
    processed = set()
    
    # Find rectangles
    for y in range(height):
        for x in range(width):
            if (y, x) in processed or grid[y][x] not in '#█':
                continue
                
            # Found a potential top-left corner
            # Look for bottom-right corner
            curr_x = x
            while curr_x < width and grid[y][curr_x] in '#█':
                curr_x += 1
            curr_x -= 1
            
            curr_y = y
            while curr_y < height and grid[curr_y][x] in '#█':
                curr_y += 1
            curr_y -= 1
            
            # Validate rectangle
            if is_valid_rect(y, x, curr_y, curr_x):
                count += 1
                # Mark the rectangle area as processed
                for py in range(y, curr_y + 1):
                    for px in range(x, curr_x + 1):
                        processed.add((py, px))
    
    return count

# Input grid
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                     ##########",
    "                                                                     #        #",
    "                     ######################################          #        #",
    "                     #                          ##########█########  #        #",
    "                     #                          #         #       #  #        #",
    "                     #                          #         #       #  #        #",
    "                     ###########################█##########       #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                        ########################█#################█##█######  #",
    "                        #                       #                 #  #     #  #",
    "                        ########################█#################█##█######  #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #      ###########█##█###     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  ###█######",
    "                                                #      #          #     #      ",
    "                                                #      #          #     #      ",
    "                                                #      #          #     #      ",
    "                                                #######█###########     #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                          #############█################█      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "       ##############                     #############█########████####█      ",
    "       #            #                                  #########█##█#####      ",
    "       ##############                                           ####           "
]

print(find_rectangles(grid))