def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        return True

    def mark_rectangle(r, c):
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and is_rectangle(x, y):
                visited[x][y] = True
                # Check all four directions
                stack.append((x + 1, y))
                stack.append((x - 1, y))
                stack.append((x, y + 1))
                stack.append((x, y - 1))

    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c):
                rectangle_count += 1
                mark_rectangle(r, c)

    return rectangle_count

# Define the grid
grid = [
    "                                               ####                             ",
    "                                               #  #                             ",
    "                                               #  #                             ",
    "                                               #  #                             ",
    "                                          #####█##█###################          ",
    "                                          #    #  #                  #          ",
    "                                          #    #  #                  #          ",
    "                                          #    #  #                  #          ",
    "                                          #    #  #                  #          ",
    "                                          #    #  #                  #          ",
    "                                          #    #  #                  #          ",
    "                                      ####█    #  #                  #          ",
    "                                      #   █    #  #                  #          ",
    "         ################███████######█## █    #  #                  #          ",
    "         #               #     #      # # █    #  #                  #          ",
    "         #               #     #      # # █    #  #                  #          ",
    "         #               #     #      # # █    #  #                  #          ",
    "         #               ######█######█## █    #  #                  #          ",
    "         #                     #      #   █    #  #                  #          ",
    "         #                     #      #   █    ####                  #          ",
    "         #                     #      #   █                          #          ",
    "         #                     #      #   █                          #          ",
    "         #                     #      #   █                          #          ",
    "         #                     #      #   █                          #          ",
    "         #                     #      #   █###########################          ",
    "         #                     #      #   #                                     ",
    "         #                     #      #   #                                     ",
    "         #                     #      #   #                                     ",
    "         #                     #      #   #                                     ",
    "         #                     #      #   #                                     ",
    "         #                     #      #   #                                     ",
    "         #                     #      #####                                     ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #                     #                                                ",
    "         #              #######█#######################                          ",
    "         #              #      #                      #                          ",
    "         #              #      #                      #                          ",
    "         #              #      #                      #                          ",
    "         ###############█#######                      #                          ",
    "                        #                             #                          ",
    "                        #           ##################█                          ",
    "                        #           #                 █                          ",
    "                        #           #                 █                          ",
    "                        #           #                 █                          ",
    "                        #           #                 █                          ",
    "                        #           #                 █                          ",
    "                        #           #                 █                          ",
    "                        ############█#################█                          ",
    "                                    #                 #                          ",
    "                                    #                 #                          ",
    "                                    #                 #                          ",
    "                                    #                 #                          ",
    "                                    #                 #                          ",
    "         ###########################█#################█##############            ",
    "         #                          #                 #             #            ",
    "         #                          #                 #             #            ",
    "         #                          #                 #             #            ",
    "         #                          #                #█##########   #            ",
    "         #                          #                ##         #   #      ####  ",
    "         #                          #                #█##########   #      #  #  ",
    "         #                          #                 #             #      #  #  ",
    "         ###########################███████████████████##############      #  #  ",
    "                                                                           ####  "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)