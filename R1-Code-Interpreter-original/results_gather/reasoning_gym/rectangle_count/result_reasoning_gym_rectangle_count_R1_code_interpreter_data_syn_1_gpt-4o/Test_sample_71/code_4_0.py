def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0
    rectangles = set()

    def find_rectangle(r, c):
        # Find the extent of the rectangle
        right = c
        while right < width and (rows[r][right] == '#' or rows[r][right] == '█'):
            right += 1
        bottom = r
        while bottom < height and (rows[bottom][c] == '#' or rows[bottom][c] == '█'):
            bottom += 1
        return (r, c, bottom - 1, right - 1)

    for i in range(height):
        for j in range(width):
            if rows[i][j] == '#' or rows[i][j] == '█':
                rect = find_rectangle(i, j)
                rectangles.add(rect)

    return len(rectangles)

grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                          ###   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                           ################################################ #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                  #####      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           ###################################█###█######## #   
                                                              #   #       # #   
                                                              #   #       # #   
                                                              #   #       █#█###
                                                              #   #       █ #  #
                                                              #   #       █ #  #
                                                              #####       █ #  #
                                                                          █ #  #
                                                                          █ #  #
                                                                          █##  #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                #################         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                       #########################█████████████████######## #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # ######
                       #                       #########                #       
                       #                       #       #                #       
                       #                       #       #                #       
             ##########█████████###############█████████#################       
             #                 #                                                
             ###################                                                
"""

print(count_rectangles(grid))