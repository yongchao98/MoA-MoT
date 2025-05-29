def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    
    def mark_rectangle(r, c):
        # Find the boundaries of the rectangle
        start_r, start_c = r, c
        end_r, end_c = r, c
        
        # Extend to the right
        while end_c + 1 < cols and grid[start_r][end_c + 1] in ['#', '█']:
            end_c += 1
        
        # Extend downwards
        while end_r + 1 < rows and all(grid[end_r + 1][col] in ['#', '█'] for col in range(start_c, end_c + 1)):
            end_r += 1
        
        # Mark all cells in the rectangle as visited
        for i in range(start_r, end_r + 1):
            for j in range(start_c, end_c + 1):
                visited[i][j] = True
    
    rectangle_count = 0
    
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ['#', '█'] and not visited[r][c]:
                mark_rectangle(r, c)
                rectangle_count += 1
    
    return rectangle_count

# Define the grid
grid = [
    "                ###################                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #                 #                                             ",
    "                #        #########█#######                                      ",
    "                #        #        #      #                                      ",
    "                #        #        #      #                                      ",
    "                #        #        #      #                                      ",
    "                #        #        #      #                                      ",
    "   #######      #        #        #      #                                      ",
    "   #     #      #        #        #      #                                      ",
    "   #     #      #        #        #      #                                      ",
    "   #     #      #        #        #      #                                      ",
    "   #     #      #        #        #      #                                      ",
    "   #     #      #        #########█#######                                      ",
    "   #     #      #                 #                                             ",
    "   #     #      #                 #                                             ",
    "   #  ###█######█#################█####################                          ",
    "   #  #  #      ###################                   #                          ",
    "   #  #  #                                       #####█###########               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                       #    #          #               ",
    "   #  #  #                                  #####█####█##########█#############  ",
    "   #  #  #                                  #    #    #          #            #  ",
    "   #  #  #                                  #    # ###█          #            #  ",
    "   #  #  #                                  #    # #  █          #            #  ",
    "   #  #  #                                  #    # #  █          #            #  ",
    "   ###█###                                  #    # #  █          #            #  ",
    "      #                                     #    # #  █      ####█############█#  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      #                                     #    # #  █      #   #            ##  ",
    "      ######################################█####█#█##█      #   #            ##  ",
    "                                            #    # #  #      #   #            ##  ",
    "                                            #    # ####      #   #            ##  ",
    "                                            #    #           #   #            ##  ",
    "                                            #    #           #   #            ##  ",
    "                           #################█####█#########  #   #            ##  ",
    "                           #                #    #        #  #   #            ##  ",
    "                           #                #    #        #  #   #            ##  ",
    "                           #                #    #        #  #   #            ##  ",
    "                           #                #    #        #  #   #            ##  ",
    "                           #                #    #        #  #   #            ##  ",
    "                   #########                #    #        #  #   #            ##  ",
    "                   #      ##                #    #        #  #   #            ##  ",
    "                   #      ##                #    #        #  #   #            ##  ",
    "                   #      ##                #    #        #  #   #            ##  ",
    "                   #      ##                #####█########█##█###█##############  ",
    "                   #      ##                     #########█##█####             #  ",
    "                   #      ##                              #  #                 #  ",
    "                   #      #################################  ###################  ",
    "                   ########                                                      "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)