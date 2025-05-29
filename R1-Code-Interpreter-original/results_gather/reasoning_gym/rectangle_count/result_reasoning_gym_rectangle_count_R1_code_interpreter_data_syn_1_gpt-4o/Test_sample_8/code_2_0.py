def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    def is_valid_rectangle(r1, c1, r2, c2):
        # Check top and bottom edges
        for c in range(c1, c2 + 1):
            if grid[r1][c] not in ('#', '█') or grid[r2][c] not in ('#', '█'):
                return False
        # Check left and right edges
        for r in range(r1, r2 + 1):
            if grid[r][c1] not in ('#', '█') or grid[r][c2] not in ('#', '█'):
                return False
        return True

    for r1 in range(rows):
        for c1 in range(cols):
            if grid[r1][c1] in ('#', '█'):
                for r2 in range(r1 + 1, rows):
                    for c2 in range(c1 + 1, cols):
                        if grid[r2][c2] in ('#', '█') and is_valid_rectangle(r1, c1, r2, c2):
                            rectangles.add((r1, c1, r2, c2))

    return len(rectangles)

# Define the grid
grid = [
    "                                                                          ###   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                                          # #   ",
    "                                                    ####                  # #   ",
    "              ###################################   #  #                  # #   ",
    "              #                                 #   #  #                  # #   ",
    "              #                                 #   #  #                  # #   ",
    "              #           ###########           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              #           #         #           #   #  #                  # #   ",
    "              ############█#########█############   #  #                  ###   ",
    "                          #         #               #  #                        ",
    "                          #         #               #  #                        ",
    "                          #         #               #  #                        ",
    "                          #         #               #  #                        ",
    "                          #         #               #  #                        ",
    "       ###################█######## #               #  #                        ",
    "       #                  #       # #               #  #                        ",
    "       #                  #       # #               #  #                        ",
    "       ###################█######## #               #  #                        ",
    "                          #         #               #  #                        ",
    "                          ###########               ####                        ",
    "                                                                 #############  ",
    "                                                        ###      #           #  ",
    "                                                        # #      #           #  ",
    "                                                        # #      #############  ",
    "                                                        # #                     ",
    "                                                        ###                     ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                      #####################################################     ",
    "                      #                                                   #     ",
    "                      #####################################################     "
]

# Count the rectangles
print(count_rectangles(grid))