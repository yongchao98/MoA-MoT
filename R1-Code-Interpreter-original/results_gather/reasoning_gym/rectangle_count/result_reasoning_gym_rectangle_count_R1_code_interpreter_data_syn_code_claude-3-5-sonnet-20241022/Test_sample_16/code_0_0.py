def is_valid_rectangle(grid, top, left, bottom, right):
    # Check horizontal edges
    for x in range(left, right + 1):
        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
            return False
    
    # Check vertical edges
    for y in range(top, bottom + 1):
        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
            return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert input to list of strings for easier processing
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            # Skip if not a corner
            if grid[top][left] not in '#█':
                continue
            
            # Try all possible bottom-right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if grid[bottom][right] not in '#█':
                        continue
                    
                    if is_valid_rectangle(grid, top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Input grid
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

print(find_rectangles(grid))