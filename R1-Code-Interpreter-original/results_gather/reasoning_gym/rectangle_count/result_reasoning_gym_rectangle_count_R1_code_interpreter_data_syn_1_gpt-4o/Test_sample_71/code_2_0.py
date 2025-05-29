def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0
    visited = [[False] * width for _ in range(height)]
    rectangle_count = 0

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        if visited[r][c]:
            return
        visited[r][c] = True
        # Find the extent of the rectangle
        right = c
        while right < width and (rows[r][right] == '#' or rows[r][right] == '█'):
            right += 1
        bottom = r
        while bottom < height and (rows[bottom][c] == '#' or rows[bottom][c] == '█'):
            bottom += 1
        # Mark the entire rectangle as visited
        for i in range(r, bottom):
            for j in range(c, right):
                visited[i][j] = True
        rectangle_count += 1

    for i in range(height):
        for j in range(width):
            if (rows[i][j] == '#' or rows[i][j] == '█') and not visited[i][j]:
                mark_rectangle(i, j)

    return rectangle_count

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