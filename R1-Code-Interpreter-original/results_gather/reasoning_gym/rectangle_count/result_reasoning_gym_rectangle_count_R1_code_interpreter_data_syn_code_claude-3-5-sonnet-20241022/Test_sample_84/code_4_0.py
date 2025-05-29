def find_rectangles(grid_str):
    # Convert string to list of lines and remove empty lines
    grid = [list(line) for line in grid_str.splitlines() if line]
    rows, cols = len(grid), len(grid[0])
    
    def is_valid_edge(y, x):
        if 0 <= y < rows and 0 <= x < cols:
            return grid[y][x] in ['#', '█']
        return False
    
    def validate_rectangle(top, left, bottom, right):
        # Skip if dimensions are too small
        if bottom <= top + 1 or right <= left + 1:
            return False
            
        # Check if all four corners exist
        if not all(is_valid_edge(y, x) for y, x in [
            (top, left), (top, right),
            (bottom, left), (bottom, right)
        ]):
            return False
            
        # Validate top and bottom edges
        for x in range(left + 1, right):
            if not is_valid_edge(top, x) or not is_valid_edge(bottom, x):
                return False
                
        # Validate left and right edges
        for y in range(top + 1, bottom):
            if not is_valid_edge(y, left) or not is_valid_edge(y, right):
                return False
                
        # Check if this is a minimal rectangle
        # (no other corners inside except at the edges)
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if (y != top and y != bottom and 
                    x != left and x != right and 
                    grid[y][x] in ['#', '█']):
                    return False
        
        return True

    # Find all corner positions
    corners = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['#', '█']:
                corners.append((i, j))
    
    rectangles = set()
    
    # Check each possible top-left corner
    for top_left in corners:
        top, left = top_left
        # Find possible bottom-right corners
        for bottom_right in corners:
            bottom, right = bottom_right
            if bottom > top and right > left:
                if validate_rectangle(top, left, bottom, right):
                    rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

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

print(find_rectangles(grid))