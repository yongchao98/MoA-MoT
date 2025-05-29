def count_rectangles(grid):
    rows = grid.strip().split('\n')
    num_rows = len(rows)
    num_cols = len(rows[0])
    
    visited = [[False] * num_cols for _ in range(num_rows)]
    rectangle_count = 0

    def is_new_rectangle(r, c):
        # Check if this is a new rectangle by looking at the surrounding
        # If the top and left are not part of a rectangle, it's a new one
        if r > 0 and (rows[r-1][c] == '#' or rows[r-1][c] == '█'):
            return False
        if c > 0 and (rows[r][c-1] == '#' or rows[r][c-1] == '█'):
            return False
        return True

    for r in range(num_rows):
        for c in range(num_cols):
            if rows[r][c] == '#' and not visited[r][c]:
                if is_new_rectangle(r, c):
                    rectangle_count += 1
                # Mark the entire rectangle as visited
                for i in range(r, num_rows):
                    if rows[i][c] != '#' and rows[i][c] != '█':
                        break
                    for j in range(c, num_cols):
                        if rows[i][j] != '#' and rows[i][j] != '█':
                            break
                        visited[i][j] = True
            elif rows[r][c] == '█' and not visited[r][c]:
                # '█' indicates overlap, count as two rectangles
                rectangle_count += 2
                # Mark the entire overlapping area as visited
                for i in range(r, num_rows):
                    if rows[i][c] != '█':
                        break
                    for j in range(c, num_cols):
                        if rows[i][j] != '█':
                            break
                        visited[i][j] = True

    return rectangle_count

# The given ASCII grid
ascii_grid = """
                       #############                                            
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
                                       #####################                    
"""

# Count the rectangles
rectangle_count = count_rectangles(ascii_grid)
print(rectangle_count)