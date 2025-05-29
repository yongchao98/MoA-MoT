def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    rectangles = 0

    def is_corner(r, c):
        if grid[r][c] == '#':
            if (r == 0 or grid[r-1][c] == ' ') and (c == 0 or grid[r][c-1] == ' '):
                return True
        return False

    def mark_rectangle(r, c):
        nonlocal rectangles
        if visited[r][c]:
            return
        visited[r][c] = True
        if grid[r][c] == '#':
            if is_corner(r, c):
                rectangles += 1
            if r + 1 < rows and grid[r+1][c] == '#':
                mark_rectangle(r+1, c)
            if c + 1 < cols and grid[r][c+1] == '#':
                mark_rectangle(r, c+1)

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#' and not visited[r][c]:
                mark_rectangle(r, c)

    return rectangles

grid = [
    "#################################################################           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                                               #           ",
    "#                                              #################█########   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                              #                #       #   ",
    "#                                   ###########████#############█########   ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             #           ",
    "#                                   #             #             ########    ",
    "####################################█#############█###############     #    ",
    "                                    #             #              #     #    ",
    "                                    #             #              #     #    ",
    "                                    #             #              #     #    ",
    "                                    #             #              #     #    ",
    "                                    #             #              #     #    ",
    "                                    #             #              #     #    ",
    "                                    ###############              #     #    ",
    "                                    #                           #######    ",
    "                                    #                            ########  ",
    "                                    #                            #      #  ",
    "                                    #                            ########  "
]

print(count_rectangles(grid))