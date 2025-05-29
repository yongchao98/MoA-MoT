def count_rectangles(grid):
    def is_corner(x, y):
        return grid[x][y] in ['#', '█']

    def trace_rectangle(x, y):
        max_x, max_y = x, y
        # Move right to find the horizontal boundary
        while max_y + 1 < len(grid[0]) and grid[x][max_y + 1] in ['#', '█']:
            max_y += 1
        # Move down to find the vertical boundary
        while max_x + 1 < len(grid) and grid[max_x + 1][y] in ['#', '█']:
            max_x += 1
        return max_x, max_y

    rectangle_count = 0
    visited = set()

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if is_corner(i, j) and (i, j) not in visited:
                bottom_right = trace_rectangle(i, j)
                # Mark all cells of the rectangle as visited
                for x in range(i, bottom_right[0] + 1):
                    for y in range(j, bottom_right[1] + 1):
                        visited.add((x, y))
                rectangle_count += 1

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