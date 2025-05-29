def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangles = set()

    def find_rectangle(r, c):
        if (r, c) in visited or grid[r][c] not in '#█':
            return None
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_r < rows and grid[end_r][c] in '#█':
            end_r += 1
        while end_c < cols and grid[r][end_c] in '#█':
            end_c += 1
        # Mark all cells in this rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited.add((i, j))
        return (r, c, end_r - 1, end_c - 1)

    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] in '#█':
                rect = find_rectangle(r, c)
                if rect:
                    rectangles.add(rect)

    return len(rectangles)

# Define the grid
grid = [
    "                                                                     ########## ",
    "                                                                     #        # ",
    "                     ######################################          #        # ",
    "                     #                          ##########█########  #        # ",
    "                     #                          #         #       #  #        # ",
    "                     #                          #         #       #  #        # ",
    "                     ###########################█##########       #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                        ########################█#################█##█######  # ",
    "                        #                       #                 #  #     #  # ",
    "                        ########################█#################█##█######  # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #      ###########█##█###     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  ###█###### ",
    "                                                #      #          #     #       ",
    "                                                #      #          #     #       ",
    "                                                #      #          #     #       ",
    "                                                #######█###########     #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                          #############█################█       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "       ##############                     #############█########████####█       ",
    "       #            #                                  #########█##█#####       ",
    "       ##############                                           ####            "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)