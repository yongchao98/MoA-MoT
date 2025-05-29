def find_rectangles(grid_str):
    # Convert string to list of lines and remove empty lines
    grid = [list(line) for line in grid_str.splitlines() if line]
    
    # Find all corner positions (# and █)
    corners = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['#', '█']:
                corners.append((i, j))
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if not all(grid[y][x] in ['#', '█'] for y, x in [(top, left), (top, right), (bottom, left), (bottom, right)]):
            return False
            
        # Check top and bottom edges
        for x in range(left + 1, right):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
                
        # Check left and right edges
        for y in range(top + 1, bottom):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
        
        return True
    
    rectangles = set()
    # Only check corners as potential top-left points
    for i, top_left in enumerate(corners):
        for j, bottom_right in enumerate(corners[i+1:], i+1):
            if bottom_right[0] > top_left[0] and bottom_right[1] > top_left[1]:  # Ensure bottom_right is actually below and right
                # Find the other two corners
                top_right = (top_left[0], bottom_right[1])
                bottom_left = (bottom_right[0], top_left[1])
                
                # Check if it forms a valid rectangle
                if is_valid_rectangle(top_left[0], top_left[1], bottom_right[0], bottom_right[1]):
                    rectangles.add((top_left[0], top_left[1], bottom_right[0], bottom_right[1]))
    
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