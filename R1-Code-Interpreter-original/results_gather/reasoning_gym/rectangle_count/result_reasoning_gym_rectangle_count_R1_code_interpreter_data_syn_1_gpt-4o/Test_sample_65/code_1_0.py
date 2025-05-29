def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0]) if num_rows > 0 else 0
    visited = [[False] * num_cols for _ in range(num_rows)]
    rectangle_count = 0

    def is_new_rectangle(r, c):
        # Check if this is a new rectangle starting point
        if rows[r][c] in '#█' and not visited[r][c]:
            # Check if it's the top-left corner of a rectangle
            if (r == 0 or rows[r-1][c] not in '#█') and (c == 0 or rows[r][c-1] not in '#█'):
                return True
        return False

    def mark_rectangle(r, c):
        # Mark the entire rectangle as visited
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < num_rows and 0 <= y < num_cols and not visited[x][y] and rows[x][y] in '#█':
                visited[x][y] = True
                stack.append((x+1, y))
                stack.append((x, y+1))

    for r in range(num_rows):
        for c in range(num_cols):
            if is_new_rectangle(r, c):
                rectangle_count += 1
                mark_rectangle(r, c)
            elif rows[r][c] == '█' and not visited[r][c]:
                # Count overlapping rectangles separately
                rectangle_count += 2
                mark_rectangle(r, c)

    return rectangle_count

# The given ASCII grid
ascii_grid = """
                                                                                
                                                                                
                                                                                
                                  ##########################                    
                                  #                        #                    
                                  #                        #           #########
                                  #                   #####█########   #       #
                                  #                   #    #       #   #       #
                                  #                   #    #       #   #       #
                                  #                   #    #       #   #########
                                  #                   #####█########            
                                  #                        #                    
                                  #                        #                    
                                  #                        #                    
                                  ######################████##############      
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
                                                        ##################      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)