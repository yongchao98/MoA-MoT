def count_rectangles(grid):
    rows = len(grid)
    cols = max(len(row) for row in grid) if rows > 0 else 0
    visited = set()
    rectangle_count = 0

    def is_rectangle_start(r, c):
        return c < len(grid[r]) and grid[r][c] in ('#', '█') and (r, c) not in visited

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < len(grid[r]) and grid[r][end_c] in ('#', '█'):
            end_c += 1
        while end_r < rows and c < len(grid[end_r]) and grid[end_r][c] in ('#', '█'):
            end_r += 1
        # Mark all positions in the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited.add((i, j))
        rectangle_count += 1

    for r in range(rows):
        for c in range(cols):
            if is_rectangle_start(r, c):
                mark_rectangle(r, c)

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