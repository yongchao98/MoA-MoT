def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangle_count = 0

    def is_rectangle_boundary(r, c):
        return grid[r][c] in ('#', '█')

    def mark_rectangle(r, c):
        if (r, c) in visited or not is_rectangle_boundary(r, c):
            return False
        # Find the extent of the rectangle
        right = c
        while right < cols and is_rectangle_boundary(r, right):
            right += 1
        down = r
        while down < rows and is_rectangle_boundary(down, c):
            down += 1
        # Mark the rectangle as visited
        for i in range(r, down):
            for j in range(c, right):
                visited.add((i, j))
        return True

    for r in range(rows):
        for c in range(cols):
            if mark_rectangle(r, c):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                                     #######    ",
    "                                                                     #   ##█    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     #   # █    ",
    "                                                                     ####█#█    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                         # #    ",
    "                                                                 ########█##    ",
    "                                                                 #       ###    ",
    "                                                                 #       ###    ",
    "#####################################################            #       ###    ",
    "#                                                   #            #       ###    ",
    "#                                                   #            #       ###    ",
    "#                                                   #            #       ###    ",
    "#                                                  #█############█#####  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                                  ##            #    #  ###    ",
    "#                                       ###########██############█####█##███    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "#                                       #          ##            #    #  ##█    ",
    "########################################█##########█#            #    #  ##█    ",
    "                                        #          #             #    #  ##█    ",
    "                                        #          #             #    #  ##█    ",
    "                                        #          #             #    #  ##█    ",
    "                                        #          #             #    #  #██    ",
    "                                        #          #             #    #  #██    ",
    "                                        #          ##############█#####  ##█    ",
    "                                        #                        #       ##█    ",
    "                                        #                        #       ##█    ",
    "                ########################█########                #       #██    ",
    "                #                       #       #                #        ##    ",
    "                #                       #       #                #        ##    ",
    "                #                       #       #                ###########    ",
    "                #                       #       #                          #    ",
    "                #                       #       #                          #    ",
    "                #                       #       #                          #    ",
    "                #                       #       #                          #    ",
    "                #                       #       #                          #    ",
    "                ########################█########                          #    ",
    "                  ######################█###########                       #    ",
    "                  #                     #          #                       #    ",
    "                  ######################█###########                       #    ",
    "                                   #####████████████████████████############    ",
    "                                   #                           #                ",
    "                                   #                           #                ",
    "                                   #                           #                ",
    "                                   #                           #                ",
    "                                   #############################                "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)