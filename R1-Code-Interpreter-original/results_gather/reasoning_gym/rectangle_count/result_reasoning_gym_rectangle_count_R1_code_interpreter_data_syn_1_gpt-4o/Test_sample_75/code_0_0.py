def count_rectangles(grid):
    def is_corner(x, y):
        # Check if the current position is a corner of a rectangle
        if grid[x][y] != '#':
            return False
        # Check surrounding positions to determine if it's a corner
        if (x > 0 and grid[x-1][y] == '#') and (y > 0 and grid[x][y-1] == '#'):
            return True
        if (x > 0 and grid[x-1][y] == '#') and (y < len(grid[0]) - 1 and grid[x][y+1] == '#'):
            return True
        if (x < len(grid) - 1 and grid[x+1][y] == '#') and (y > 0 and grid[x][y-1] == '#'):
            return True
        if (x < len(grid) - 1 and grid[x+1][y] == '#') and (y < len(grid[0]) - 1 and grid[x][y+1] == '#'):
            return True
        return False

    def mark_rectangle(x, y):
        # Mark the rectangle starting from (x, y)
        stack = [(x, y)]
        while stack:
            cx, cy = stack.pop()
            if grid[cx][cy] == '#':
                grid[cx][cy] = '.'
                if cx > 0 and grid[cx-1][cy] == '#':
                    stack.append((cx-1, cy))
                if cx < len(grid) - 1 and grid[cx+1][cy] == '#':
                    stack.append((cx+1, cy))
                if cy > 0 and grid[cx][cy-1] == '#':
                    stack.append((cx, cy-1))
                if cy < len(grid[0]) - 1 and grid[cx][cy+1] == '#':
                    stack.append((cx, cy+1))

    rectangle_count = 0
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if is_corner(i, j):
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

# Define the grid
grid = [
    "                                   #######################                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #                     #                      ",
    "                                   #######################                      ",
    "                                                                          ######",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                                                                          #    #",
    "                 ##############################                           #    #",
    "                 #                            #                           #    #",
    "                 #                            #                           #    #",
    "                 ##############################                           ######",
    "                                #####################                           ",
    "                                #                   #                           ",
    "                                #####################                           "
]

# Convert grid to list of lists
grid = [list(row) for row in grid]

# Count rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)