def count_rectangles(grid):
    rows = [line.rstrip() for line in grid.strip().split('\n')]
    height = len(rows)
    width = max(len(row) for row in rows) if height > 0 else 0

    visited = [[False] * width for _ in range(height)]

    def is_new_rectangle(x, y):
        if visited[x][y]:
            return False
        if rows[x][y] == '#':
            if (x == 0 or (y < len(rows[x-1]) and rows[x-1][y] != '#')) and (y == 0 or rows[x][y-1] != '#'):
                return True
        return False

    def mark_rectangle(x, y):
        stack = [(x, y)]
        while stack:
            cx, cy = stack.pop()
            if 0 <= cx < height and 0 <= cy < len(rows[cx]) and not visited[cx][cy]:
                if rows[cx][cy] in '#█':
                    visited[cx][cy] = True
                    if cx + 1 < height:
                        stack.append((cx+1, cy))
                    if cy + 1 < len(rows[cx]):
                        stack.append((cx, cy+1))
                    if cx - 1 >= 0:
                        stack.append((cx-1, cy))
                    if cy - 1 >= 0:
                        stack.append((cx, cy-1))

    rectangle_count = 0

    for i in range(height):
        for j in range(len(rows[i])):
            if is_new_rectangle(i, j):
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

grid = """
   ################                                                             
   #              #                                                             
   #              #                                                             
   #              #                                                             
   #              #                                                             
   #              #                              ###############                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #                
   #              #                              #             #         #####  
   #              #                              #             #         #   #  
   #              #                              #             #         #   #  
   ################                              #             #         #   #  
                                                 #             #         #   #  
                                                 #             #         #   #  
                                                 #             #         #   #  
                                                 #             #         #   #  
                                                 #           ##█         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █         #   #  
                                                 #           # █ #########   #  
                                                 #           # █ #      ##   #  
                                                 #           # █ #      ##   #  
                                                 #           # █ #      ##   #  
                                                 ############█#█ #      ##   #  
                                                             # # #      ##   #  
                                                             # # #      ##   #  
                                                             # # #      ##   #  
                                                             # # #      ##   #  
                                                             # # #########   #  
                                                             # #         #####  
                                                             # #                
                                                             # #                
                                                             # #                
        #####################################################█#█##              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    # # #              
        #                                                    ### #              
        ##########################################################              
"""

print(count_rectangles(grid))