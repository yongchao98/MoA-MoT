def count_rectangles(grid):
    def is_rectangle(x1, y1, x2, y2):
        # Check if the boundary forms a rectangle
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
        return True

    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangle_count = 0
    visited = set()

    # Iterate over all possible top-left corners of rectangles
    for y1 in range(rows):
        for x1 in range(cols):
            if grid[y1][x1] in '#█' and (x1, y1) not in visited:
                # Find the bottom-right corner
                for y2 in range(y1, rows):
                    for x2 in range(x1, cols):
                        if grid[y2][x2] in '#█' and is_rectangle(x1, y1, x2, y2):
                            # Mark all positions in this rectangle as visited
                            for y in range(y1, y2 + 1):
                                for x in range(x1, x2 + 1):
                                    visited.add((x, y))
                            rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                               #####            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                         ###############       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #####            ",
    "                                         #             #       ## ##            ",
    "                                         #             #       ## ##            ",
    "                                         #             #       #####            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         #             #       #   #            ",
    "                                         ###############       #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #   #            ",
    "                                                               #####            ",
    "                                                                                 ",
    "              ########################################                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #                                      #                          ",
    "              #               ############           #                          ",
    "              #               #          #           #                          ",
    "              #               #          #           #                          ",
    "              #               #          #           #                          ",
    "              #               #          #           #                          ",
    "              #               #          #           #                          ",
    "              #               #          #           #                          ",
    "              #               ############           #                          ",
    "              ########################################                          "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)