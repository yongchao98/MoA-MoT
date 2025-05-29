def is_rectangle(grid, x1, y1, x2, y2):
    # Check if the region forms a valid rectangle
    for x in range(x1, x2 + 1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
    for y in range(y1, y2 + 1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Convert grid to list of lists for easier processing
    grid = [list(row) for row in grid]
    
    # Find all potential top-left corners
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] not in '#█':
                continue
                
            # For each top-left corner, try to find bottom-right corner
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if grid[y2][x2] not in '#█':
                        continue
                        
                    # Check if this forms a valid rectangle
                    if is_rectangle(grid, x1, y1, x2, y2):
                        # Check if it's a minimal rectangle (no smaller rectangle inside)
                        is_minimal = True
                        for y in range(y1 + 1, y2):
                            for x in range(x1 + 1, x2):
                                if grid[y][x] in '#█':
                                    is_minimal = False
                                    break
                            if not is_minimal:
                                break
                        if is_minimal:
                            rectangles.append((x1, y1, x2, y2))

    return len(rectangles)

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
    "                                                            #############       ",
    "                                          ##################█############       ",
    "                   #######################█#                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      ##                #          ##       ",
    "                   #                      #█################█############       ",
    "                   #                       #                #           #       ",
    "                   #                       #                #           #       ",
    "                   #                    ###█#############   #           #       ",
    "                   #                    #  #            #   #           #       ",
    "                   #                    #  #            #   #           #       ",
    "                   #                    #  #            #   #           #       ",
    "                   #                    #  #            #   #           #       ",
    "                   #                    #  #            #   #           #       ",
    "                   #                    ###█#############   #           #       ",
    "                   #              #########█#################           #       ",
    "                   #              #        #               ##           #       ",
    "                   #              #        #               ##           #       ",
    "                   #         #####█########█#########      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##           #       ",
    "                   #         #    #        #        #      ##############       ",
    "                   #         #    #        #        #      #                    ",
    "                   #         #    #        #        #      #                    ",
    "                   #         #    #        #        #      #                 ###",
    "                   #         #    #        #        #      #                 # #",
    "                   #         #    #        #        #      #                 ###",
    "                   ##########█####█#########        #      #                    ",
    "                             #    #                 #      #                    ",
    "                             #    #                 #      #                    ",
    "                             #    #                 #      #                    ",
    "                             #    #                 #      #                    ",
    "                             #####█##################      #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                            ######█########################█##############      ",
    "                            #     #                        #             #      ",
    "                            ######█########################█##############      ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  #                        #                    ",
    "                                  ##########################                    ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))