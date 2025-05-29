def count_rectangles(grid):
    rows = grid.strip().split('\n')
    height = len(rows)
    width = max(len(row) for row in rows)
    
    visited = set()
    rectangles = 0

    def is_new_rectangle(x, y):
        if (x, y) in visited:
            return False
        if x >= len(rows[y]) or rows[y][x] not in '#█':
            return False
        # Check if it's a top-left corner of a rectangle
        if x > 0 and x < len(rows[y]) and rows[y][x-1] in '#█':
            return False
        if y > 0 and x < len(rows[y-1]) and rows[y-1][x] in '#█':
            return False
        return True

    for y in range(height):
        for x in range(len(rows[y])):
            if is_new_rectangle(x, y):
                rectangles += 1
                # Mark the whole rectangle as visited
                x_end = x
                while x_end < len(rows[y]) and rows[y][x_end] in '#█':
                    x_end += 1
                y_end = y
                while y_end < height and x < len(rows[y_end]) and rows[y_end][x] in '#█':
                    y_end += 1
                for i in range(x, x_end):
                    for j in range(y, y_end):
                        visited.add((i, j))

    return rectangles

grid = """
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

print(count_rectangles(grid))