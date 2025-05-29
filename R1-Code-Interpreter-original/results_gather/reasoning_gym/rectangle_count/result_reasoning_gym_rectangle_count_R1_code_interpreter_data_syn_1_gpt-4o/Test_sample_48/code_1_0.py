# Define the grid as a list of strings
grid = [
    "                                                                   ############ ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   #          # ",
    "                                                                   ############ ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                      ######### ",
    "                                                                      #       # ",
    "      ####################################################            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      #                                                  #            #       # ",
    "      ####################################################            #       # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                              ########################█###    # ",
    "                                              #                       #  #    # ",
    "                                              ########################█###    # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                                  #####################       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #                  ##       # ",
    "                                                  #####################       # ",
    "                                                                      #       # ",
    "                                                                      #       # ",
    "                                                                      ######### ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

# Function to count rectangles
def count_rectangles(grid):
    rectangles = 0
    visited = set()
    
    def is_corner(x, y):
        return grid[x][y] == '#' or grid[x][y] == '█'
    
    def find_rectangle(x, y):
        if (x, y) in visited or not is_corner(x, y):
            return None
        visited.add((x, y))
        # Find the bottom-right corner
        max_x, max_y = x, y
        while max_x + 1 < len(grid) and is_corner(max_x + 1, y):
            max_x += 1
        while max_y + 1 < len(grid[0]) and is_corner(x, max_y + 1):
            max_y += 1
        # Mark all corners of the rectangle
        for i in range(x, max_x + 1):
            for j in range(y, max_y + 1):
                visited.add((i, j))
        return (x, y, max_x, max_y)
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if is_corner(i, j) and (i, j) not in visited:
                if find_rectangle(i, j):
                    rectangles += 1
    
    return rectangles

# Count the rectangles
total_rectangles = count_rectangles(grid)
print(total_rectangles)