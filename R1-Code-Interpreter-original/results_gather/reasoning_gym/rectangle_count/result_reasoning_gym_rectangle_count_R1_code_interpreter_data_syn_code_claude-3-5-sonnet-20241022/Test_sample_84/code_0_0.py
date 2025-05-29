def is_corner(grid, i, j):
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in ['#', '█']

def is_edge(grid, i, j):
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in ['#', '█']

def verify_rectangle(grid, top, left, bottom, right):
    # Check corners
    if not all(is_corner(grid, i, j) for i, j in [
        (top, left), (top, right), (bottom, left), (bottom, right)
    ]):
        return False
    
    # Check horizontal edges
    for j in range(left + 1, right):
        if not (is_edge(grid, top, j) and is_edge(grid, bottom, j)):
            return False
    
    # Check vertical edges
    for i in range(top + 1, bottom):
        if not (is_edge(grid, i, left) and is_edge(grid, i, right)):
            return False
    
    return True

def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows, cols = len(grid), len(grid[0])
    rectangles = 0
    
    # Convert input to 2D array
    grid = [list(row) for row in grid.split('\n')]
    
    # Try all possible rectangle combinations
    for top in range(rows):
        for left in range(cols):
            for bottom in range(top + 1, rows):
                for right in range(left + 1, cols):
                    if verify_rectangle(grid, top, left, bottom, right):
                        rectangles += 1
    
    return rectangles

# Input grid
grid = """                                                                                
                                                                     #######    
                                                                     #   ##█    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     #   # █    
                                                                     ####█#█    
                                                                         # #    
                                                                         # #    
                                                                         # #    
                                                                         # #    
                                                                         # #    
                                                                         # #    
                                                                         # #    
                                                                 ########█##    
                                                                 #       ###    
                                                                 #       ###    
#####################################################            #       ###    
#                                                   #            #       ###    
#                                                   #            #       ###    
#                                                   #            #       ###    
#                                                  #█############█#####  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                                  ##            #    #  ###    
#                                       ###########██############█####█##███    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
#                                       #          ##            #    #  ##█    
########################################█##########█#            #    #  ##█    
                                        #          #             #    #  ##█    
                                        #          #             #    #  ##█    
                                        #          #             #    #  ##█    
                                        #          #             #    #  #██    
                                        #          #             #    #  #██    
                                        #          ##############█#####  ##█    
                                        #                        #       ##█    
                                        #                        #       ##█    
                ########################█########                #       #██    
                #                       #       #                #        ##    
                #                       #       #                #        ##    
                #                       #       #                ###########    
                #                       #       #                          #    
                #                       #       #                          #    
                #                       #       #                          #    
                #                       #       #                          #    
                #                       #       #                          #    
                ########################█########                          #    
                  ######################█###########                       #    
                  #                     #          #                       #    
                  ######################█###########                       #    
                                   #####████████████████████████############    
                                   #                           #                
                                   #                           #                
                                   #                           #                
                                   #                           #                
                                   #############################                """

print(count_rectangles(grid))