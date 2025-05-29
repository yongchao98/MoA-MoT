def find_rectangles(grid_str):
    # Convert string to list of lines and remove empty lines
    grid = [list(line) for line in grid_str.splitlines() if line]
    rows, cols = len(grid), len(grid[0])
    
    def is_complete_rectangle(top, left, bottom, right):
        # Check if this forms a valid rectangle by verifying:
        # 1. All corners are present
        # 2. All edges are continuous
        # 3. No other corners exist within the edges
        
        # Check corners
        corners = [(top, left), (top, right), (bottom, left), (bottom, right)]
        if not all(grid[y][x] in ['#', '█'] for y, x in corners):
            return False
        
        # Check horizontal edges and ensure no extra corners in between
        for x in range(left + 1, right):
            # Top edge
            if grid[top][x] not in ['#', '█']:
                return False
            # Bottom edge
            if grid[bottom][x] not in ['#', '█']:
                return False
        
        # Check vertical edges and ensure no extra corners in between
        for y in range(top + 1, bottom):
            # Left edge
            if grid[y][left] not in ['#', '█']:
                return False
            # Right edge
            if grid[y][right] not in ['#', '█']:
                return False
        
        return True

    def find_next_corner(start_row, start_col):
        # Find the next '#' or '█' from the given position
        for i in range(start_row, rows):
            start = start_col if i == start_row else 0
            for j in range(start, cols):
                if grid[i][j] in ['#', '█']:
                    return (i, j)
        return None

    rectangles = set()
    
    # Find top-left corners
    for top in range(rows):
        for left in range(cols):
            if grid[top][left] not in ['#', '█']:
                continue
                
            # Find possible right edges from this left edge
            for right in range(left + 1, cols):
                if grid[top][right] not in ['#', '█']:
                    continue
                    
                # Find possible bottom edges
                for bottom in range(top + 1, rows):
                    if grid[bottom][left] not in ['#', '█'] or grid[bottom][right] not in ['#', '█']:
                        continue
                        
                    # Verify if this forms a valid rectangle
                    if is_complete_rectangle(top, left, bottom, right):
                        # Check if this rectangle is minimal (no smaller rectangles within)
                        is_minimal = True
                        for y in range(top + 1, bottom):
                            for x in range(left + 1, right):
                                if grid[y][x] in ['#', '█']:
                                    is_minimal = False
                                    break
                            if not is_minimal:
                                break
                        
                        if is_minimal:
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