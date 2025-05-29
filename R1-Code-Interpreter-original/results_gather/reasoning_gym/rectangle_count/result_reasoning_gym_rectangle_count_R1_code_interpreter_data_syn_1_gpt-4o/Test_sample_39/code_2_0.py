def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0

    visited = [[False] * width for _ in range(height)]
    rectangle_count = 0

    def is_rectangle_start(x, y):
        return (rows[x][y] == '#' or rows[x][y] == '█') and not visited[x][y]

    def mark_rectangle(x, y):
        nonlocal rectangle_count
        if visited[x][y]:
            return
        visited[x][y] = True
        if rows[x][y] == '█':
            rectangle_count += 1
        elif rows[x][y] == '#':
            # Trace the rectangle
            right = y
            while right < width and (rows[x][right] == '#' or rows[x][right] == '█'):
                right += 1
            bottom = x
            while bottom < height and (rows[bottom][y] == '#' or rows[bottom][y] == '█'):
                bottom += 1
            # Mark the rectangle
            for i in range(x, bottom):
                for j in range(y, right):
                    if rows[i][j] == '#' or rows[i][j] == '█':
                        visited[i][j] = True
            rectangle_count += 1

    for i in range(height):
        for j in range(width):
            if is_rectangle_start(i, j):
                mark_rectangle(i, j)

    return rectangle_count

grid = """
                                                                 ############   
                                                                 #          #   
                                                                 #          #   
                                            #############        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                                            #           #        #          #   
                   #########################█###########█########█########  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        #       #  #   
                   #                        #           #        ########█###   
                   #                        #           #                #      
                   #                        #           #                #      
                   #                        #           #                #      
                   #                        #           #                #      
                   #                        #           #                #      
                   #                        #           #                #      
                   #########################█###########█#################      
                                            #           #                       
                                            #           #                       
                                            #           #                       
                                            #           #                       
                                            #           #                       
                                            #           #                       
                                            #############                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                   ########     
                                                                   #      #     
                                                                   #      #     
                                                                   #      #     
                                                                   #      #     
                                                                   #      #     
                                                                   ######## ### 
           ############################                                     # # 
           #                          #                                     ### 
           ############################                                         
"""

print(count_rectangles(grid))