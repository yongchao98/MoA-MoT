def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def mark_rectangle(r, c):
        # Mark the entire rectangle as visited
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and not visited[x][y] and grid[x][y] in ('#', '█'):
                visited[x][y] = True
                # Check all 4 directions
                stack.append((x + 1, y))
                stack.append((x - 1, y))
                stack.append((x, y + 1))
                stack.append((x, y - 1))

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