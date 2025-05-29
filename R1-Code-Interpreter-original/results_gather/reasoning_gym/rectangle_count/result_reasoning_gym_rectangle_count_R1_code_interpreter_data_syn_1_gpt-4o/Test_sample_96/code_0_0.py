def count_rectangles(grid):
    def is_corner(x, y):
        # Check if the current position is a corner of a rectangle
        if grid[x][y] not in ['#', '█']:
            return False
        # Check surrounding positions to determine if it's a corner
        if (x > 0 and grid[x-1][y] in [' ', '#', '█']) and (y > 0 and grid[x][y-1] in [' ', '#', '█']):
            return True
        if (x < len(grid) - 1 and grid[x+1][y] in [' ', '#', '█']) and (y < len(grid[0]) - 1 and grid[x][y+1] in [' ', '#', '█']):
            return True
        return False

    def find_rectangle(x, y):
        # Find the bottom-right corner of the rectangle starting from (x, y)
        max_x, max_y = x, y
        while max_x < len(grid) and grid[max_x][y] in ['#', '█']:
            max_x += 1
        while max_y < len(grid[0]) and grid[x][max_y] in ['#', '█']:
            max_y += 1
        return max_x - 1, max_y - 1

    visited = set()
    rectangle_count = 0

    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and is_corner(i, j):
                br_x, br_y = find_rectangle(i, j)
                # Mark all positions of this rectangle as visited
                for x in range(i, br_x + 1):
                    for y in range(j, br_y + 1):
                        visited.add((x, y))
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                      #######                                   ",
    "                                      #     #                                   ",
    "                                      #     #                                   ",
    "                                      #     #                                   ",
    "                                      #     #                                   ",
    "                                      #######                ###                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             ###                ",
    "                                                   ###################          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   #                 #          ",
    "                                                   ###################          ",
    "                                                                     #######    ",
    "                                                                     #     #    ",
    "                                                                     #     #    ",
    "                                                                     #     #    ",
    "                                                                     #     #    ",
    "                                                                     #     #    ",
    "                                                                     #######    "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)