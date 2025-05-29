def count_rectangles(grid):
    rows = len(grid)
    cols = max(len(row) for row in grid) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    
    def is_rectangle(r, c):
        if r >= rows or c >= len(grid[r]):
            return False
        if grid[r][c] not in ['#', '█'] or visited[r][c]:
            return False
        return True
    
    def mark_rectangle(r, c):
        # Find the bottom-right corner of the rectangle
        start_r, start_c = r, c
        while c < len(grid[start_r]) and grid[r][c] in ['#', '█']:
            c += 1
        c -= 1
        while r < rows and c < len(grid[r]) and grid[r][c] in ['#', '█']:
            r += 1
        r -= 1
        
        # Mark the rectangle as visited
        for i in range(start_r, r + 1):
            for j in range(start_c, c + 1):
                visited[i][j] = True
    
    rectangle_count = 0
    for r in range(rows):
        for c in range(len(grid[r])):
            if is_rectangle(r, c):
                mark_rectangle(r, c)
                rectangle_count += 1
    
    return rectangle_count

# Define the grid
grid = [
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "  ###########                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                          ############################             ",
    "  #         #                          #                          #             ",
    "  #         #                          #                          #             ",
    "  #         #                          ############################             ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                                                                    ",
    "  #         #                    #########################                      ",
    "  #         #                    #                       #                      ",
    "  #         #                    #                       #                      ",
    "  #         #                    #                       #                      ",
    "  #         #       #############█############           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  #         #       #            #           #           #                      ",
    "  ###########       #            #           #           #                      ",
    "                    #            #           #           #                      ",
    "                    #            ############█############                      ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #                                  ",
    "                    #                        #            ###                   ",
    "                    #                        #            # #                   ",
    "                    #                       #█############█#█############       ",
    "                    #                       ##            # #           #       ",
    "                    #                       ##            # #           #       ",
    "                    #                       ##            # #           #       ",
    "                    #                       ##            # #           #       ",
    "                    #                       ##            # #           #       ",
    "                    #                       ##            # #           #       ",
    "                    #         ##############██############█#█###########█###### ",
    "                    #         #             ##            # #           # ####█ ",
    "                    #         #             ##            # #           # #   █ ",
    "                    #         #             #█############█#█############ #   █ ",
    "                    #         #              #            # #             #   █ ",
    "                    ##########█###############            # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           # #             #   █ ",
    "                              #                           ###             #   █ ",
    "                              #                                           #   █ ",
    "                              #                                           #   █ ",
    "                              #                                           #   █ ",
    "                              ############################################█████ ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)