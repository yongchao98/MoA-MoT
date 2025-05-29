def find_rectangles(grid):
    # Split the grid into lines
    lines = grid.strip().split('\n')
    height = len(lines)
    width = len(lines[0])
    
    def is_corner(y, x):
        # Check if position is a corner of a rectangle
        if y < 0 or y >= height or x < 0 or x >= width:
            return False
        if lines[y][x] not in '#█':
            return False
        
        # Count adjacent '#' or '█' characters
        adjacent = 0
        if y > 0 and lines[y-1][x] in '#█':
            adjacent += 1
        if y < height-1 and lines[y+1][x] in '#█':
            adjacent += 1
        if x > 0 and lines[y][x-1] in '#█':
            adjacent += 1
        if x < width-1 and lines[y][x+1] in '#█':
            adjacent += 1
        
        return adjacent == 2

    def find_rectangle(start_y, start_x):
        # Find the bottom-right corner of the rectangle
        right_x = start_x
        while right_x < width and lines[start_y][right_x] in '#█':
            right_x += 1
        right_x -= 1
        
        bottom_y = start_y
        while bottom_y < height and lines[bottom_y][start_x] in '#█':
            bottom_y += 1
        bottom_y -= 1
        
        # Verify it's a valid rectangle
        for y in range(start_y, bottom_y + 1):
            for x in range(start_x, right_x + 1):
                if y in (start_y, bottom_y) or x in (start_x, right_x):
                    if lines[y][x] not in '#█':
                        return None
                else:
                    if lines[y][x] != ' ':
                        return None
        
        return (start_y, start_x, bottom_y, right_x)

    rectangles = []
    # Find all rectangles by looking for top-left corners
    for y in range(height):
        for x in range(width):
            if (lines[y][x] in '#█' and 
                (y == 0 or lines[y-1][x] not in '#█') and 
                (x == 0 or lines[y][x-1] not in '#█')):
                rect = find_rectangle(y, x)
                if rect:
                    rectangles.append(rect)

    # Count overlapping points
    overlap_count = sum(1 for y in range(height) 
                       for x in range(width) 
                       if lines[y][x] == '█')

    total_rectangles = len(rectangles) + overlap_count

    print(total_rectangles)

# Input grid
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