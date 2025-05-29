def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_corner(r, c):
        if grid[r][c] not in ('#', '█'):
            return False
        if r > 0 and grid[r-1][c] in ('#', '█'):
            return False
        if c > 0 and grid[r][c-1] in ('#', '█'):
            return False
        return True

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        if visited[r][c]:
            return
        visited[r][c] = True
        if grid[r][c] in ('#', '█'):
            if is_corner(r, c):
                rectangle_count += 1
            if r + 1 < rows and grid[r+1][c] in ('#', '█'):
                mark_rectangle(r+1, c)
            if c + 1 < cols and grid[r][c+1] in ('#', '█'):
                mark_rectangle(r, c+1)

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ('#', '█') and not visited[r][c]:
                mark_rectangle(r, c)

    return rectangle_count

# Define the grid
grid = [
    "                    ####################################                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #                        ",
    "                    #                                  #    ####                ",
    "         ###########█##################################█####█##█#######         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    #  #      #         ",
    "         #          #                                  #    ####      #         ",
    "         #          #                                  #              #         ",
    "         #          #                                  #              #         ",
    "         #     #####█#####                             #              #         ",
    "         #     #    #    #                             #              #         ",
    "         #     #    #    #                             #              #         ",
    "         ######█####█####█#############################█###############         ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #                        ",
    "               #    #    #                             #              ########  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              #      #  ",
    "               #    #    #                             #              ########  ",
    "               #    #    #          ###################█#                       ",
    "               #    #    #          #                  ##                       ",
    "               #    #    #          #                  #█################       ",
    "               #    #    #          #                  #█               #       ",
    "               #####█#####          #                  #█               #       ",
    "                    ################█████████████████████               #       ",
    "                                                        #               #       ",
    "                                                        #               #       ",
    "                                                        #               #       ",
    "                                                        #               #       ",
    "                                                        #               #       ",
    "                                                        #################       "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)