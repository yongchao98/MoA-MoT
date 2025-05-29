def find_rectangles(grid):
    # Convert grid to 2D array
    grid = [list(row) for row in grid.splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_corner(i, j):
        if grid[i][j] not in ['#', '█']:
            return False
        # Count adjacent '#' or '█' characters
        neighbors = 0
        if i > 0 and grid[i-1][j] in ['#', '█']: neighbors += 1
        if i < height-1 and grid[i+1][j] in ['#', '█']: neighbors += 1
        if j > 0 and grid[i][j-1] in ['#', '█']: neighbors += 1
        if j < width-1 and grid[i][j+1] in ['#', '█']: neighbors += 1
        return neighbors == 2

    def trace_rectangle(start_i, start_j):
        # Find the opposite corner
        end_i = start_i
        end_j = start_j
        
        # Find right edge
        j = start_j
        while j < width and grid[start_i][j] in ['#', '█']:
            j += 1
        right_j = j - 1
        
        # Find bottom edge
        i = start_i
        while i < height and grid[i][start_j] in ['#', '█']:
            i += 1
        bottom_i = i - 1
        
        # Verify rectangle
        if not is_corner(bottom_i, right_j):
            return False
            
        # Mark visited
        for i in range(start_i, bottom_i + 1):
            for j in range(start_j, right_j + 1):
                if grid[i][j] in ['#', '█']:
                    grid[i][j] = 'V'
        return True

    rectangle_count = 0
    # Scan for top-left corners
    for i in range(height):
        for j in range(width):
            if grid[i][j] in ['#', '█'] and is_corner(i, j):
                if trace_rectangle(i, j):
                    rectangle_count += 1

    return rectangle_count

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                           #############################################        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                               ########    #        
                           #                               #      #    #        
                       ####█#####                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          ########    #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
         ##############█###█####█##########################         ####        
         #             ####█#####                         #         # ##        
         #                 #                              #         ####        
         #                 #                           ###█##########  #        
         #                 ############################█##█#########█###        
         # #########################################   #  #         #           
         # #                                       #   #  #         #           
         # #                                       #   #  #         #           
         # #                                       #   █###         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #########################################   █###         #           
         #                                             ###█##########           
         #                                                #                     
         #                                                #                     
         #                                                #                     
         #                                                #                     
         #                             ###################█##########           
         ##############################█###################         #           
                                       #                            #           
                                       ##############################           """

print(find_rectangles(grid))