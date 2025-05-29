def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    visited = []

    # Initialize visited with the correct width for each row
    for row in rows:
        visited.append([False] * len(row))

    rectangle_count = 0

    def is_new_rectangle(x, y):
        if visited[x][y]:
            return False
        if rows[x][y] == '#':
            # Check if it's a top-left corner of a rectangle
            if (x == 0 or (y < len(rows[x-1]) and rows[x-1][y] != '#')) and (y == 0 or rows[x][y-1] != '#'):
                return True
        elif rows[x][y] == '█':
            # Check if it's part of overlapping rectangles
            if (x == 0 or (y < len(rows[x-1]) and rows[x-1][y] != '█')) and (y == 0 or rows[x][y-1] != '█'):
                return True
        return False

    for i in range(height):
        width = len(rows[i])
        for j in range(width):
            if is_new_rectangle(i, j):
                rectangle_count += 1
                # Mark the entire rectangle as visited
                stack = [(i, j)]
                while stack:
                    x, y = stack.pop()
                    if 0 <= x < height and 0 <= y < len(rows[x]) and not visited[x][y] and rows[x][y] in '#█':
                        visited[x][y] = True
                        # Add neighbors to the stack
                        if x + 1 < height and y < len(rows[x + 1]):
                            stack.append((x + 1, y))
                        if y + 1 < len(rows[x]):
                            stack.append((x, y + 1))
                        if x - 1 >= 0 and y < len(rows[x - 1]):
                            stack.append((x - 1, y))
                        if y - 1 >= 0:
                            stack.append((x, y - 1))

    return rectangle_count

grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                    ####                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    #  #                        
                                                    ####                        
                                                        ###                     
                                                        # #                     
                                                        # #                     
                                                        # #                     
                                                        # #                     
                                                        # #                     
                                                        ###                     
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                    ########                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
                                                    #      #                    
      ##############################################█######█##############      
      #                                             #      #             #      
      #                                             #      #             #      
      #                                             #      #             #      
      #                                             #      #             #      
      #                                             #      #             #      
      #                                             #      #             #      
      #                                             #      #             #      
######█                                             #      #             #      
#     █                                             #      #             #      
#     █                                             #      #             #      
#     █                                             #      #             #      
#     █                                             #      #             #      
#     █#############################################████████##############      
#     #                                                                         
#     #                                                                         
#######                                                                         
                                                                                
"""

print(count_rectangles(grid))