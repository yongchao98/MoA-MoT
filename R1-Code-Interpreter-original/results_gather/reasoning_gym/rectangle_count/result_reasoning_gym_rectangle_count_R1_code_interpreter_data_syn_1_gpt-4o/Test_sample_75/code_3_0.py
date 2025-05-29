def count_rectangles(grid):
    rectangles = 0
    visited = set()
    
    def is_corner(x, y):
        return (grid[y][x] == '#' or grid[y][x] == '█') and \
               (x == 0 or grid[y][x-1] == ' ') and \
               (y == 0 or grid[y-1][x] == ' ')
    
    def mark_rectangle(x, y):
        width = 0
        height = 0
        # Find width
        while x + width < len(grid[0]) and (grid[y][x + width] == '#' or grid[y][x + width] == '█'):
            width += 1
        # Find height
        while y + height < len(grid) and (grid[y + height][x] == '#' or grid[y + height][x] == '█'):
            height += 1
        # Mark all positions of this rectangle as visited
        for i in range(height):
            for j in range(width):
                visited.add((x + j, y + i))
    
    for y in range(len(grid)):
        for x in range(len(grid[0])):
            if is_corner(x, y) and (x, y) not in visited:
                rectangles += 1
                mark_rectangle(x, y)
    
    return rectangles

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

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)