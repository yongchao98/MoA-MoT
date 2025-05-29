def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangles = []

    def is_top_left_corner(r, c):
        # Check if this is a top-left corner of a rectangle
        if grid[r][c] in ('#', '█'):
            if (r == 0 or grid[r-1][c] == ' ') and (c == 0 or grid[r][c-1] == ' '):
                return True
        return False

    def find_rectangle(r, c):
        # Find the bottom-right corner of the rectangle starting at (r, c)
        if not is_top_left_corner(r, c):
            return None
        # Find the right boundary
        right = c
        while right < cols and grid[r][right] in ('#', '█'):
            right += 1
        # Find the bottom boundary
        bottom = r
        while bottom < rows and grid[bottom][c] in ('#', '█'):
            bottom += 1
        # Check if this forms a valid rectangle
        for i in range(r, bottom):
            if grid[i][right-1] not in ('#', '█'):
                return None
        for j in range(c, right):
            if grid[bottom-1][j] not in ('#', '█'):
                return None
        return (r, c, bottom-1, right-1)

    for r in range(rows):
        for c in range(cols):
            if is_top_left_corner(r, c):
                rect = find_rectangle(r, c)
                if rect and rect not in rectangles:
                    rectangles.append(rect)

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