def count_rectangles(grid):
    # Function to check if a cell is a corner of a rectangle
    def is_corner(x, y):
        return grid[x][y] in ['#', '█'] and (
            (grid[x-1][y] in ['#', '█'] and grid[x][y-1] in ['#', '█']) or
            (grid[x+1][y] in ['#', '█'] and grid[x][y+1] in ['#', '█']) or
            (grid[x-1][y] in ['#', '█'] and grid[x][y+1] in ['#', '█']) or
            (grid[x+1][y] in ['#', '█'] and grid[x][y-1] in ['#', '█'])
        )

    # Parse the grid to find all corners
    corners = []
    for i in range(1, len(grid) - 1):
        for j in range(1, len(grid[i]) - 1):
            if is_corner(i, j):
                corners.append((i, j))

    # Count rectangles by finding pairs of corners
    rectangle_count = 0
    visited = set()
    for (x1, y1) in corners:
        for (x2, y2) in corners:
            if x1 < x2 and y1 < y2:
                if (x1, y2) in corners and (x2, y1) in corners:
                    if (x1, y1, x2, y2) not in visited:
                        rectangle_count += 1
                        visited.add((x1, y1, x2, y2))
                        visited.add((x2, y2, x1, y1))

    return rectangle_count

# Define the grid
grid = [
    "               ############################################################     ",
    "               #                                                          #     ",
    "               ############################################################     ",
    "                        ###########████████████████████████#####################",
    "                        #          #                      #                    #",
    "                        #          #                      #                    #",
    "                        ###########█######################█#####################",
    "                                   #                      #                     ",
    "                                   #                      #                     ",
    "                                   #                      #                     ",
    "                                   ########################                     "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)