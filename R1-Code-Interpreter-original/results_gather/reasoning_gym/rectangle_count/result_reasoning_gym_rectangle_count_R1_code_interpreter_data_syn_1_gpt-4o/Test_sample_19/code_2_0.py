def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def mark_rectangle(r, c):
        # Find the bottom-right corner of the rectangle
        max_r, max_c = r, c
        while max_r < rows and grid[max_r][c] in ('#', '█'):
            max_r += 1
        while max_c < cols and grid[r][max_c] in ('#', '█'):
            max_c += 1
        # Mark the entire rectangle as visited
        for i in range(r, max_r):
            for j in range(c, max_c):
                if grid[i][j] in ('#', '█'):
                    visited[i][j] = True

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ('#', '█') and not visited[i][j]:
                # Found a new rectangle
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

# Define the grid
grid = [
    "                                                                ###########     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                                                                #         #     ",
    "                             ###################################█##       #     ",
    "                             #                                  # #       #     ",
    "                             #                                  # #       #     ",
    "                             #                                  # #       #     ",
    "                             #                                  # #       #     ",
    "                             #                                  ##█########     ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #             ",
    "                             #                                    #        #####",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             #                                    #        #   #",
    "                             ######################################        #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                           #   #",
    "                                                                      #####█###█",
    "                                                                      #    #   █",
    "                                                                      #    ####█",
    "                                                                      #        #",
    "                                                                      #        #",
    "                                                                      #        #",
    "                                                                      #        #",
    "                                                                      ##########"
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)