def is_minimal_rectangle(grid, top, left, bottom, right):
    # Check if this is a minimal rectangle (not containing other rectangles)
    
    # First verify it's a valid rectangle outline
    # Check horizontal edges
    for x in range(left, right + 1):
        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
            return False
    
    # Check vertical edges
    for y in range(top, bottom + 1):
        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
            return False
    
    # Check that interior is mostly empty (except for overlaps)
    for y in range(top + 1, bottom):
        for x in range(left + 1, right):
            if grid[y][x] not in ' █':
                return False
    
    # Check that there's no rectangle boundary inside
    has_internal_boundary = False
    for y in range(top + 1, bottom):
        for x in range(left + 1, right):
            if grid[y][x] == '#':
                has_internal_boundary = True
                break
    
    return not has_internal_boundary

def find_minimal_rectangles(grid):
    # Split the grid into lines
    lines = [line for line in grid.split('\n') if line.strip()]
    grid = [list(row) for row in lines]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = set()
    
    # Find starting points (top-left corners)
    for top in range(height - 1):
        for left in range(width - 1):
            if grid[top][left] not in '#█':
                continue
                
            # Find matching bottom-right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if grid[bottom][right] not in '#█':
                        continue
                    
                    # Check if this forms a minimal rectangle
                    if is_minimal_rectangle(grid, top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
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

print(find_minimal_rectangles(grid))