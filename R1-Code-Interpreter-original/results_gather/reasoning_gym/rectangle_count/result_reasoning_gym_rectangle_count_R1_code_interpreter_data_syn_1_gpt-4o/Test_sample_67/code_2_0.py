def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = len(rows[0]) if height > 0 else 0

    visited = [[False] * width for _ in range(height)]
    rectangle_count = 0

    def is_new_rectangle(x, y):
        if visited[x][y]:
            return False
        if rows[x][y] in '#█':
            if (x == 0 or rows[x-1][y] not in '#█') and (y == 0 or rows[x][y-1] not in '#█'):
                return True
        return False

    def mark_rectangle(x, y):
        stack = [(x, y)]
        while stack:
            cx, cy = stack.pop()
            if 0 <= cx < height and 0 <= cy < width and not visited[cx][cy] and rows[cx][cy] in '#█':
                visited[cx][cy] = True
                stack.append((cx+1, cy))
                stack.append((cx, cy+1))
                stack.append((cx-1, cy))
                stack.append((cx, cy-1))

    for i in range(height):
        for j in range(width):
            if is_new_rectangle(i, j):
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                             ##################                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
           ##################█##############################################    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #               ##█#########      ###                           #    
           #               # #        #      ###                           #    
           #               # #        #      ###                           #    
           #               # #        #      ###                    ####   #    
           #               # #        #      ###                    #  #   #    
           #               # #        #      ###                    #  #   #    
           #               # #        #      ###                    #  #   #    
           #               ##█#########      ###                    #  #   #    
           ##################█████████████████##                    #  #   #    
                                               #                    #  #   #    
                                               #                    #  #   #    
                                               #                    ####   #    
                                               #                           #    
                                               #                           #    
                                               #                           #    
                                               #                           #    
                                               #     ######################█### 
                                               ######█###                  #  # 
                                               ##    ###█##################█### 
                                               #█████████###################    
                                                                                
"""

print(count_rectangles(grid))