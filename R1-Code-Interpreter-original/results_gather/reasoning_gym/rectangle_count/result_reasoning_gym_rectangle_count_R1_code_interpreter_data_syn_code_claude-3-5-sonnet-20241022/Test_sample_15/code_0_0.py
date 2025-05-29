def is_corner(grid, i, j):
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in '#█'

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid]
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if not is_corner(grid, top, left):
                continue
                
            # Find possible bottom and right coordinates
            for bottom in range(top + 1, height):
                if not is_corner(grid, bottom, left):
                    continue
                    
                for right in range(left + 1, width):
                    if not (is_corner(grid, top, right) and is_corner(grid, bottom, right)):
                        continue
                        
                    # Verify if it forms a valid rectangle
                    valid = True
                    # Check horizontal edges
                    for x in range(left + 1, right):
                        if not (grid[top][x] in '#█' and grid[bottom][x] in '#█'):
                            valid = False
                            break
                    
                    # Check vertical edges
                    for y in range(top + 1, bottom):
                        if not (grid[y][left] in '#█' and grid[y][right] in '#█'):
                            valid = False
                            break
                    
                    if valid:
                        rectangles.add((top, left, bottom, right))
    
    # Count additional rectangles from overlapping areas
    extra_rectangles = 0
    for i in range(height):
        for j in range(width):
            if grid[i][j] == '█':
                extra_rectangles += 1
                
    return len(rectangles) + extra_rectangles

# Input grid
grid = [
    "   ################                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                              ###############                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #         #####  ",
    "   #              #                              #             #         #   #  ",
    "   #              #                              #             #         #   #  ",
    "   ################                              #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #           ##█         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █ #########   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 ############█#█ #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #########   #  ",
    "                                                             # #         #####  ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "        #####################################################█#█##              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    ### #              ",
    "        ##########################################################              "
]

result = find_rectangles(grid)
print(result)