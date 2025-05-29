def count_rectangles(grid):
    def is_rectangle(x, y):
        # Check if the current position is a top-left corner of a rectangle
        if grid[x][y] not in ('#', '█'):
            return False
        # Check for a rectangle by moving right and down
        width, height = 0, 0
        while y + width < len(grid[0]) and grid[x][y + width] in ('#', '█'):
            width += 1
        while x + height < len(grid) and grid[x + height][y] in ('#', '█'):
            height += 1
        # Verify the rectangle by checking the bottom and right edges
        for i in range(width):
            if grid[x + height - 1][y + i] not in ('#', '█'):
                return False
        for i in range(height):
            if grid[x + i][y + width - 1] not in ('#', '█'):
                return False
        return True

    count = 0
    visited = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and is_rectangle(i, j):
                count += 1
                # Mark the rectangle as visited
                width, height = 0, 0
                while j + width < len(grid[0]) and grid[i][j + width] in ('#', '█'):
                    width += 1
                while i + height < len(grid) and grid[i + height][j] in ('#', '█'):
                    height += 1
                for x in range(i, i + height):
                    for y in range(j, j + width):
                        visited.add((x, y))
    return count

# Define the grid
grid = [
    "                                                     #####                      ",
    "                                                     #   #                      ",
    "                                                     #####                      ",
    "           #################################################################### ",
    "           #                                                                  # ",
    "           #                                                                  # ",
    "           #                                                                  # ",
    "           #                                                     ########     # ",
    "           #                                                     #      #     # ",
    "           #                                                     #      #     # ",
    "           #                                                     #      #     # ",
    "           #                                                     #      #     # ",
    "           #                                            #########█######█##   # ",
    "           #                                            #        #      # #   # ",
    "           #                                            #        #      # #   # ",
    "           #                                            #        #      # #   # ",
    "           #                                            #        #      # #   # ",
    "           #                                          ##█########█##### # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #    # # #   # ",
    "           #                                          # #        #####█## #   # ",
    "           #                                          # #             #   #   # ",
    "           #                                          # #             #   #   # ",
    "           #                                          # #             #   #   # ",
    "           #                                          # #             #   #   # ",
    "           #                                          # #             #   #   # ",
    "           #                                          # #             #   #   # ",
    "   ########█###################################       # #             #   #   # ",
    "   #       #                                  #       # #            #█###██##█#",
    "   #       #                                  #       # #            ##   #█  ##",
    "   #       #                                  #       # #            ##   #█##██",
    "   #       #                                  #       # #            ##   #█# #█",
    "   #       #                                  #       # #            ##   #██#██",
    "   #       #                                  #       # #            ##   ### ##",
    "   #  #####█########                          #       # #            ##   ### ##",
    "   #  #    ########█##########################█#######█#█############██###███###",
    "   #  #            #                          #       # #            ##   ###  #",
    "   #  #            #                          #       # #            ##   ###  #",
    "   #  ##############                          #       # #            ##   ###  #",
    "   ############################################       # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      # #            ##   ###  #",
    "                                                      ##█############█#   ###  #",
    "                                                        #            #    ###  #",
    "                                                        #############█##########",
    "                                                                     #     #    ",
    "                                                                     #######    ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)