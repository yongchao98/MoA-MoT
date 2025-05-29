def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangles = set()

    def is_corner(r, c):
        return grid[r][c] in ('#', '█')

    def find_rectangle(r, c):
        if (r, c) in visited or not is_corner(r, c):
            return None
        visited.add((r, c))
        # Find the bottom-right corner of the rectangle
        br_r, br_c = r, c
        while br_r + 1 < rows and grid[br_r + 1][c] in ('#', '█'):
            br_r += 1
        while br_c + 1 < cols and grid[r][br_c + 1] in ('#', '█'):
            br_c += 1
        # Mark all cells in the rectangle as visited
        for i in range(r, br_r + 1):
            for j in range(c, br_c + 1):
                visited.add((i, j))
        return (r, c, br_r, br_c)

    for r in range(rows):
        for c in range(cols):
            if is_corner(r, c) and (r, c) not in visited:
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
    "                           ################################################ #   ",
    "                           #                                             ## #   ",
    "                           #                                             ## #   ",
    "                           #                                             ## #   ",
    "                           #                                             ## #   ",
    "                           #                                             ## #   ",
    "                           #                                             ## #   ",
    "                           #                                  #####      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           #                                  #   #      ## #   ",
    "                           ###################################█###█######## #   ",
    "                                                              #   #       # #   ",
    "                                                              #   #       # #   ",
    "                                                              #   #       █#█###",
    "                                                              #   #       █ #  #",
    "                                                              #   #       █ #  #",
    "                                                              #####       █ #  #",
    "                                                                          █ #  #",
    "                                                                          █ #  #",
    "                                                                          █##  #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                #################         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                                                #               #         #    #",
    "                       #########################█████████████████######## #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # #    #",
    "                       #                                                # ######",
    "                       #                       #########                #       ",
    "                       #                       #       #                #       ",
    "                       #                       #       #                #       ",
    "             ##########█████████###############█████████#################       ",
    "             #                 #                                                 ",
    "             ###################                                                 "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)