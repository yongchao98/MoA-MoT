def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    visited = []

    # Initialize the visited matrix with the correct dimensions for each row
    for row in rows:
        visited.append([False] * len(row))

    rectangle_count = 0

    def is_new_rectangle(r, c):
        # Check if the current position is the top-left corner of a new rectangle
        if rows[r][c] == '#' and not visited[r][c]:
            if (r == 0 or c >= len(rows[r-1]) or rows[r-1][c] != '#') and (c == 0 or rows[r][c-1] != '#'):
                return True
        return False

    def mark_rectangle(r, c):
        # Mark the entire rectangle as visited
        stack = [(r, c)]
        while stack:
            x, y = stack.pop()
            if 0 <= x < num_rows and 0 <= y < len(rows[x]) and not visited[x][y] and rows[x][y] in '#█':
                visited[x][y] = True
                # Add adjacent cells to the stack
                if x + 1 < num_rows and y < len(rows[x+1]) and rows[x+1][y] in '#█' and not visited[x+1][y]:
                    stack.append((x+1, y))
                if y + 1 < len(rows[x]) and rows[x][y+1] in '#█' and not visited[x][y+1]:
                    stack.append((x, y+1))

    for r in range(num_rows):
        for c in range(len(rows[r])):
            if is_new_rectangle(r, c):
                rectangle_count += 1
                mark_rectangle(r, c)

    return rectangle_count

# ASCII grid as a string
ascii_grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                      ######################                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #                                    
                      #                    #              #####                 
                      #                    #              #   #                 
                      ######################              #   #                 
                                                          #   #                 
                                                          #   #                 
                             #############################█#  #                 
                             #                            ##  #                 
                             #                   #########██##█#######          
                             #                   #        ##  #      #          
                             ####################█########█#  #      # #########
                                                 #        #   #      # #       #
    #########                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #       #                                    #        #   #      # #       #
    #########                                    #        #   #      # #       #
                                                 #        #####      # #       #
                                                 #                   # #       #
                                                 #                   # #       #
                                                 #  #######          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #       #
                                                 #  #     #          # #########
                                                 #  #     #          #          
   #####                                         #  #     #          #          
   #   #                                         ###█#####█###########          
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #     #                     
   #   #                                            #######                     
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #   #                                                                        
   #####                                                                        
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)