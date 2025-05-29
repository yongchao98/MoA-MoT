def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = 0
    
    # Helper function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Helper function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if not is_border(grid[top][left]):
                continue
            
            # Find possible bottom and right coordinates
            for bottom in range(top + 1, height):
                if not is_border(grid[bottom][left]):
                    continue
                for right in range(left + 1, width):
                    if not is_border(grid[bottom][right]) or not is_border(grid[top][right]):
                        continue
                    
                    # Validate the rectangle
                    if validate_rectangle(top, left, bottom, right):
                        rectangles += 1

    return rectangles

# Input grid
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                       ######                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                               ##############      ####█####█###############    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               #            #      #   #    #              #    ",
    "                               ##############      #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   #   #    #              #    ",
    "                                                   ####█####█###############    ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       #    #                   ",
    "                                                       ######                   ",
    "                                                                                ",
    "                       ##################################################       ",
    "                       #                                                #       ",
    "                       #                                                #       ",
    "                       #                                                #       ",
    "                       #                                                #       ",
    "                       #                                                #       ",
    "                       #  #################################             #       ",
    "                       #  #                               #             #       ",
    "                       #  #################################             #       ",
    "                       #                                                #       ",
    "                       ##################################################       ",
    "                                                                                "
]

print(find_rectangles(grid))