def count_rectangles(grid):
    def is_top_left(x, y):
        # Check if the current position is a top-left corner of a rectangle
        if grid[x][y] not in ('#', '█'):
            return False
        # Ensure it's the top-left by checking above and left
        if x > 0 and grid[x-1][y] in ('#', '█'):
            return False
        if y > 0 and grid[x][y-1] in ('#', '█'):
            return False
        return True

    def find_rectangle(x, y):
        # Find the extent of the rectangle starting from (x, y)
        width, height = 0, 0
        while y + width < len(grid[0]) and grid[x][y + width] in ('#', '█'):
            width += 1
        while x + height < len(grid) and grid[x + height][y] in ('#', '█'):
            height += 1
        return width, height

    count = 0
    visited = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and is_top_left(i, j):
                width, height = find_rectangle(i, j)
                count += 1
                # Mark the rectangle as visited
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