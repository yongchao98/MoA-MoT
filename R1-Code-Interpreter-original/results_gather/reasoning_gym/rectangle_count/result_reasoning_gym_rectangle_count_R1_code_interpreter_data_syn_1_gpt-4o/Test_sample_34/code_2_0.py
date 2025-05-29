def count_rectangles(grid):
    rectangles = 0
    visited = set()
    
    for i, line in enumerate(grid):
        for j, char in enumerate(line):
            if char == '#' and (i, j) not in visited:
                # Found a new rectangle, let's explore it
                rectangles += 1
                # Mark all parts of this rectangle as visited
                # Find the bottom-right corner
                end_i, end_j = i, j
                while end_j < len(line) and grid[i][end_j] == '#':
                    end_j += 1
                while end_i < len(grid) and all(grid[end_i][k] == '#' for k in range(j, end_j)):
                    end_i += 1
                # Mark the rectangle as visited
                for x in range(i, end_i):
                    for y in range(j, end_j):
                        visited.add((x, y))
    
    return rectangles

# The grid as a list of strings
grid = [
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                       ##############################           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       #                            #           ",
    "                                       ##############################           ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                     ###################                        ",
    "                                     #                 #                        ",
    "                                     ###################                        ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
    "                                      ",
]

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)