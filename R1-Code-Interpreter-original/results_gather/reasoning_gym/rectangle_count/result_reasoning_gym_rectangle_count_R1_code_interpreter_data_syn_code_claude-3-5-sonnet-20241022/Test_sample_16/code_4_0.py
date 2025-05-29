def find_complete_rectangle(grid, start_y, start_x):
    height = len(grid)
    width = len(grid[0])
    
    # Find right edge
    right_x = start_x
    while right_x < width and grid[start_y][right_x] in '#█':
        right_x += 1
    right_x -= 1
    
    # Find bottom edge
    bottom_y = start_y
    while bottom_y < height and grid[bottom_y][start_x] in '#█':
        bottom_y += 1
    bottom_y -= 1
    
    # Verify rectangle
    # Check if all corners exist
    if not (grid[start_y][right_x] in '#█' and 
            grid[bottom_y][start_x] in '#█' and 
            grid[bottom_y][right_x] in '#█'):
        return None
    
    # Check all edges are complete
    for x in range(start_x, right_x + 1):
        if grid[start_y][x] not in '#█' or grid[bottom_y][x] not in '#█':
            return None
    
    for y in range(start_y, bottom_y + 1):
        if grid[y][start_x] not in '#█' or grid[y][right_x] not in '#█':
            return None
    
    # Check interior is valid (space or overlap)
    for y in range(start_y + 1, bottom_y):
        for x in range(start_x + 1, right_x):
            if grid[y][x] not in ' █':
                return None
    
    return (start_y, start_x, bottom_y, right_x)

def count_rectangles(grid):
    # Convert grid to list of lists
    grid = [list(line) for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    rectangles = set()
    
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if this could be a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    rect = find_complete_rectangle(grid, y, x)
                    if rect:
                        rectangles.add(rect)
    
    return len(rectangles)

# Input grid (same as before)
grid = """                                                                                
                                                                                
                                                                                
                                                            #########           
                                                            #       #           
                                                            #       #           
                                                            #       #           
                    ###########                             #       #           
                    #         #                             #       #           
                    #         #                             #       #           
                    #         #                             #       #           
                    #         #                             #       #           
                    #         #                             #       #           
                    #         #                             #       #           
                    #         #                             #########           
                    #         #                                                 
                    #         #                       #####                     
                    #         #                       #   #                     
                    #         #                       #   #                     
                    #         #                       ####█####                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
                    #         #                       ##  #   #                 
             #######█#########█#######################██##█###█################ 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
             #      #         #                       ##  #   #               # 
 ############█######█#########█#######################██##█###█############## # 
 #           #      #         #                       ##  #   #             # # 
 #           #      #         #                       ##  #   #             # # 
 #           #      #         #                       ##  #   #             # # 
 #           #      #         #                       ##  #   #             # # 
 #           #      #         #                       ##  #   #             # # 
 #           #      #         #                       ##  #   #             # # 
 #           #######█#########█#######################██##█###█#############█## 
 #                  #         #                       ##  #   #             #   
 #                  #         #                       ##  #   #             #   
 #                  #         #                       ##  #   #             #   
 #                  #         #                       ##  #   #             #   
 #                  #         #                       ##  #   #             #   
 #                  #         #                       #█###   #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #             #   
 #                  #         #                        #      #     ####### #   
 #                  #         #                        #      #  ###█#####█##   
 #                  #         #                        #      #  #  #     ###   
 #                  #         #                        #      #  #  #     ###   
 #                  #         #                        #      #  #  #     ###   
 #                  #         #                        #      #  ###███████##   
 #                  #         #                        #      #             #   
 #                  #         #         ###########    #      #             #   
 #                  #         #         #         #    ########             #   
 #                  #         #         #         #                         #   
 #                  #         #         #         #                         #   
 #                  #         #         #         #                         #   
 #                  #         #         #         #                         #   
 #                  #         #         #         #                         #   
 ###################█#########█#########███████████##########################   
                    ###########                                                 """

print(count_rectangles(grid))