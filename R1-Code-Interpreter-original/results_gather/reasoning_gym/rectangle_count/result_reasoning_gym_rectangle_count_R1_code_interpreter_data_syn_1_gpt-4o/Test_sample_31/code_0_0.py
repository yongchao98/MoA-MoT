def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    
    def is_rectangle(r, c):
        if grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        return True
    
    def mark_rectangle(r, c):
        # Find the extent of the rectangle
        start_r, start_c = r, c
        end_r, end_c = r, c
        
        # Find the bottom-right corner of the rectangle
        while end_c < cols and grid[start_r][end_c] in ('#', '█'):
            end_c += 1
        end_c -= 1
        
        while end_r < rows and grid[end_r][start_c] in ('#', '█'):
            end_r += 1
        end_r -= 1
        
        # Mark the rectangle as visited
        for i in range(start_r, end_r + 1):
            for j in range(start_c, end_c + 1):
                visited[i][j] = True
        
        # Count the rectangle
        return 1 if grid[start_r][start_c] == '#' else 2
    
    rectangle_count = 0
    
    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c):
                rectangle_count += mark_rectangle(r, c)
    
    return rectangle_count

# Define the grid
grid = [
    "                         #####################                                  ",
    "           ##############█###################█#################                 ",
    "           #             #                   #                #                 ",
    "           #             #                   #                #                 ",
    "           #             #                   #                #                 ",
    "           #             #                   #                #                 ",
    "           #             #                   #                #                 ",
    "       ####█#############█############       #                #                 ",
    "       #   #             #           #       #                #                 ",
    "       #   #             #           #       #                #                 ",
    "       #   #             #           #       #                #                 ",
    "       #   #             #           #       #                #                 ",
    "       #   ##############█###########█#######█#################                 ",
    "       #                 #           #       #                                  ",
    "       #                 #           #       #                                  ",
    "       #                 #           #       #                                  ",
    "       #                 #           #       #                                  ",
    "       #                 #           #       ###################################",
    "       #                 #           #       ##                                #",
    "       #                 #           #       ##                                #",
    "    ###█#################█#          #       ##                                #",
    "    #  #                 ##          #       ##                                #",
    "    #  #                 ##          #       ##                                #",
    "    #  #                 ##          #       ##         ############           #",
    "    #  #                 ##          #       ##         #          #           #",
    "    #  #                 ##          #       ###########█##########█############",
    "    #  #                 ##          #       #          #          #            ",
    "    #  #                 ##          #       #          #          #            ",
    "    #  #                 ##          #       #          #          #            ",
    "    #  #                 ##          #       #          #          #            ",
    "    #  ##################██###########       #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    ##                  #          #          #            ",
    "    #                    #█###################          #          #            ",
    "    #                     #                             #          #            ",
    "    #                     #                             #          #            ",
    "    #                     #                             #          #            ",
    "    #                     #                             #          #            ",
    "    #                     #                             #          #   ######   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             #          #   #    #   ",
    "    #                     #                             ############   #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #                                            #    #   ",
    "    #                     #   ####################                     #    #   ",
    "    #                     #   #                  #                     ######   ",
    "    #                     #   ####################                              ",
    "    #######################                                                      "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)