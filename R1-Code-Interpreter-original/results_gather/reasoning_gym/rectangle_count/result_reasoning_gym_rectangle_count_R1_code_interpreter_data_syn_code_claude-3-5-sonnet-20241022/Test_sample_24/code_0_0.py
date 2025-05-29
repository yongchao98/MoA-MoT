def is_corner(grid, i, j):
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in '#█'

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all potential top-left corners
    for i in range(height):
        for j in range(width):
            if not is_corner(grid, i, j):
                continue
                
            # Look for bottom-right corners
            for y in range(i + 1, height):
                for x in range(j + 1, width):
                    if not is_corner(grid, y, x):
                        continue
                        
                    # Check if top-right and bottom-left corners exist
                    if not (is_corner(grid, i, x) and is_corner(grid, y, j)):
                        continue
                    
                    # Verify the lines connecting corners
                    valid = True
                    # Check horizontal lines
                    for col in range(j + 1, x):
                        if not (grid[i][col] in '#█' and grid[y][col] in '#█'):
                            valid = False
                            break
                    
                    # Check vertical lines
                    for row in range(i + 1, y):
                        if not (grid[row][j] in '#█' and grid[row][x] in '#█'):
                            valid = False
                            break
                    
                    if valid:
                        rectangles.add((i, j, y, x))
    
    return len(rectangles)

# Input grid
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                 ###########################                                    ",
    "                 #                         #                                    ",
    "                 #                         #                                    ",
    "                 #                         #                             ####   ",
    "                 #                         #                             #  #   ",
    "                 #                         #   ######################    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 ###########################   #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                ###############█####################█##  #  #   ",
    "                                #              #                    # #  #  #   ",
    "                                #              #                    # #  #  #   ",
    "                                ###############█####################███##█##█#  ",
    "                                               #                    ##   #  ##  ",
    "                                               # ###################██###█##██##",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  #####█##█# #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # ###################█####█##█###",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #        ########    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #########█######█#####    #  #   ",
    "                                                        #      #         #  #   ",
    "                                                        ########         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    " ###############################################################         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #          ############################                       #         #  #   ",
    " ###########█##########################█########################         #  #   ",
    "            #                          #                                 #  #   ",
    "            ############################                                 #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         ####   "
]

print(find_rectangles(grid))