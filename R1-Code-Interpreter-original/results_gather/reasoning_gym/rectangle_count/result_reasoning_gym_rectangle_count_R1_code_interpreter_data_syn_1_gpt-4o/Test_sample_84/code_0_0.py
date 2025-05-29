def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle_boundary(r, c):
        return grid[r][c] in ('#', '█')

    def mark_rectangle(r, c):
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and not visited[x][y] and is_rectangle_boundary(x, y):
                visited[x][y] = True
                # Check all 4 directions
                stack.append((x + 1, y))
                stack.append((x - 1, y))
                stack.append((x, y + 1))
                stack.append((x, y - 1))

    for r in range(rows):
        for c in range(cols):
            if is_rectangle_boundary(r, c) and not visited[r][c]:
                mark_rectangle(r, c)
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