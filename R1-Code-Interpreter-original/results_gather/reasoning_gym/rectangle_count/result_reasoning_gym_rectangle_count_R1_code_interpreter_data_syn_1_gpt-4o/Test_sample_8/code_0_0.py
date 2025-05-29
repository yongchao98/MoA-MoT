def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    def find_rectangle(r, c):
        if grid[r][c] not in ('#', '█'):
            return None
        # Find the bottom-right corner of the rectangle
        br_r, br_c = r, c
        while br_r + 1 < rows and grid[br_r + 1][c] in ('#', '█'):
            br_r += 1
        while br_c + 1 < cols and grid[r][br_c + 1] in ('#', '█'):
            br_c += 1
        # Check if it's a valid rectangle
        for i in range(r, br_r + 1):
            if grid[i][c] not in ('#', '█') or grid[i][br_c] not in ('#', '█'):
                return None
        for j in range(c, br_c + 1):
            if grid[r][j] not in ('#', '█') or grid[br_r][j] not in ('#', '█'):
                return None
        return (r, c, br_r, br_c)

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ('#', '█'):
                rect = find_rectangle(r, c)
                if rect:
                    rectangles.add(rect)

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