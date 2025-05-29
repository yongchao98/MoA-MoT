def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle_start(r, c):
        return (grid[r][c] == '#' or grid[r][c] == '█') and not visited[r][c]

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        # Find the extent of the rectangle
        end_r, end_c = r, c
        while end_c < cols and (grid[r][end_c] == '#' or grid[r][end_c] == '█'):
            end_c += 1
        while end_r < rows and (grid[end_r][c] == '#' or grid[end_r][c] == '█'):
            end_r += 1
        # Mark the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited[i][j] = True
        rectangle_count += 1

    for r in range(rows):
        for c in range(cols):
            if is_rectangle_start(r, c):
                mark_rectangle(r, c)

    return rectangle_count

# Define the grid
grid = [
    "                                                               ###              ",
    "                                                               # #              ",
    "                                                               # #              ",
    "                                                               # #              ",
    "                                                               # #              ",
    "                                                               # #              ",
    "                                                               # #              ",
    "                     #######################################   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #              ",
    "                     #                                     #   # #       ###    ",
    "                     #                                     #   # #       # #    ",
    "                     #                                     #   # #       # #    ",
    "                     #                                     #   # #       # #    ",
    "                     #                   ########          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #                   #      #          #   # #       # #    ",
    "                     #              #####█######█##########█###█#█       # #    ",
    "                     #              #    #      #          #   # █       # #    ",
    "                     #              #    #      #          #   ##█       # #    ",
    "       ##############█##############█####█######█##########█#####█#######█#█#   ",
    "       #             #              #    #      #          #     #       # ##   ",
    "       #             #              #    #      #          #     #       # ##   ",
    "       #             #              #    #      #          #     #       # ##   ",
    "       #             #              #    #      #          #     #       # ##   ",
    "       #             #              #    #      #          #     #       # ##   ",
    "       #             #              #    #      #       ###█#####█#######█ ##   ",
    "       #             #              #    #      #       #  #     #       █ ##   ",
    "       #             #              #    #███████####   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #   ######█## #  #     #       █ ##   ",
    "       #             #              #     #   #     # # #  #     #       █ ##   ",
    "       #             #              #     #   ######█## #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     #         #   #  #     #       █ ##   ",
    "       #             #              #     ###########   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   #  #     #       #  #   ",
    "       #             #              #                   ###█#####█########  #   ",
    "       #             #              #                      #     #          #   ",
    "       #             #              #######################█######          #   ",
    "       #             #                                     #                #   ",
    "       #             #######################################                #   ",
    "       #                                                                    #   ",
    "       #                                                                    #   ",
    "       #                                                                    #   ",
    "       #                                                                    #   ",
    "       #                                                                    #   ",
    "       ######################################################################   "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)