def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangles = []

    def is_new_rectangle(r, c):
        # Check if this is a new rectangle starting point
        if grid[r][c] in ('#', '█') and not visited[r][c]:
            # Check if it's a top-left corner of a rectangle
            if (r == 0 or grid[r-1][c] == ' ') and (c == 0 or grid[r][c-1] == ' '):
                return True
        return False

    def mark_rectangle(r, c):
        # Mark the rectangle as visited
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and not visited[x][y] and grid[x][y] in ('#', '█'):
                visited[x][y] = True
                # Add adjacent cells
                stack.append((x+1, y))
                stack.append((x, y+1))
                stack.append((x-1, y))
                stack.append((x, y-1))

    for r in range(rows):
        for c in range(cols):
            if is_new_rectangle(r, c):
                rectangles.append((r, c))
                mark_rectangle(r, c)

    return len(rectangles)

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