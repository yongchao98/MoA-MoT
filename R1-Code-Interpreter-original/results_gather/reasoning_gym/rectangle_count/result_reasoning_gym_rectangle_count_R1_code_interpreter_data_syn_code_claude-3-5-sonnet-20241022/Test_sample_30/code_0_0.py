def find_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid.strip().split('\n')]
    height = len(grid)
    width = max(len(row) for row in grid)
    
    # Pad grid to make all rows equal length
    for i in range(height):
        grid[i].extend([' '] * (width - len(grid[i])))

    def is_corner(i, j):
        return grid[i][j] in '#█'

    def trace_rectangle(start_i, start_j):
        # Find right edge
        right_j = start_j
        while right_j < width and grid[start_i][right_j] in '#█':
            right_j += 1
        right_j -= 1

        # Find bottom edge
        bottom_i = start_i
        while bottom_i < height and grid[bottom_i][start_j] in '#█':
            bottom_i += 1
        bottom_i -= 1

        # Verify rectangle
        if not all(grid[start_i][j] in '#█' for j in range(start_j, right_j + 1)):
            return False
        if not all(grid[bottom_i][j] in '#█' for j in range(start_j, right_j + 1)):
            return False
        if not all(grid[i][start_j] in '#█' for i in range(start_i, bottom_i + 1)):
            return False
        if not all(grid[i][right_j] in '#█' for i in range(start_i, bottom_i + 1)):
            return False

        return True

    rectangles = 0
    visited = set()

    # Scan for top-left corners
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):
                # Check if this could be a top-left corner
                if (i == 0 or grid[i-1][j] not in '#█') and (j == 0 or grid[i][j-1] not in '#█'):
                    if trace_rectangle(i, j):
                        rectangles += 1

    print(rectangles)

# Test with the provided grid
grid = """                       #############                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                     ##█###########█########################################### 
                     # #           #                                          # 
                     # #           #                                          # 
                     ##█###########█########################################### 
                       #           #                                            
                       #           #                                            
                       #           #   #########################                
                       #           #   #                       #                
                       #           #   #                       #                
                       #           #   #                       #                
                       #           #   #                       #                
                       #           #   #                       #                
                       #           #   #   #############       #                
                       #           #   #   #           #       #                
                       #           #   #   #           #       #                
                      #█#####      #   #   #           #       #                
                      ##    #      #   #   #           #       #                
                      ##    #      #   #   #           #       #                
                      ##    #      #   #   #           #       #                
                      ##    #      #   #   #           #       #                
                      ######█#######   #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #                
                      #     #          #   #           #       #        ###     
                      #     #          #   #           #       #        # #     
                      #     #          #   #           #       #        # #     
                      #     #          #   #           #       #        # #     
                      #     #          #   #           #       #        # #     
                      #     #          ####█###########█########        # #     
                      #     #              #           #                # #     
                      #     #              #           #                # #     
                      #     #              #           #                # #     
                      #     #              █###########█#############   # #     
                      #     #              █           #            #   # #     
                      #######              █###########█#############   ###     
                                           #           #                        
                                           #           #                        
                                           #           #                        
                                           #           #                        
                                           #           #    ##################  
                                           #           #    #                #  
                                           #           #    ##################  
          #################################█########## #                        
          #                                #         # #                        
          #                                #         # #                        
          #                                #         # #                        
          #                                #         # #                        
          #                                ##########█##                        
          #                                          #                          
          #                                          #                          
          #                            ##############█######                    
          #############################█##############     #                    
                                       #####################                    """

find_rectangles(grid)