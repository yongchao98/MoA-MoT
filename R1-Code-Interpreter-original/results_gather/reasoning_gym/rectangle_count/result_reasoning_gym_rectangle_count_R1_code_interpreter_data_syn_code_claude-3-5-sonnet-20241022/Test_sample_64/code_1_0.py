def is_minimal_rectangle(grid, top, left, bottom, right):
    # Check if this is a minimal rectangle (no subdivisions inside)
    
    # Check corners
    if grid[top][left] not in '#█' or grid[top][right] not in '#█' or \
       grid[bottom][left] not in '#█' or grid[bottom][right] not in '#█':
        return False
    
    # Check horizontal edges
    for x in range(left + 1, right):
        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
            return False
        
    # Check vertical edges
    for y in range(top + 1, bottom):
        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
            return False
    
    # Check inside - should be empty except for overlapping points
    for y in range(top + 1, bottom):
        for x in range(left + 1, right):
            if grid[y][x] not in ' █':
                return False
    
    # Check if this rectangle is part of a larger one
    has_extension = False
    
    # Check if can extend left
    if left > 0 and all(grid[y][left-1] in '#█' for y in range(top, bottom + 1)):
        has_extension = True
    # Check if can extend right
    if right < len(grid[0])-1 and all(grid[y][right+1] in '#█' for y in range(top, bottom + 1)):
        has_extension = True
    # Check if can extend up
    if top > 0 and all(grid[top-1][x] in '#█' for x in range(left, right + 1)):
        has_extension = True
    # Check if can extend down
    if bottom < len(grid)-1 and all(grid[bottom+1][x] in '#█' for x in range(left, right + 1)):
        has_extension = True
        
    return not has_extension

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in '#█':
                continue
            
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if is_minimal_rectangle(grid, top, left, bottom, right):
                        count += 1
    
    return count

# Process the input grid
grid = """
                ###################################                             
                #                          #####  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #   #  #                             
                #                          #####  #                             
                #                       ##########█#################            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       #         #                #            
                #                       ##########█#################            
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 # #######################     
                #                                 # #                     #     
                #                                 # #                     #     
                #                                 # #                     #     
                #                                 # #                     #     
                #               ##################█#█###                  #     
                #               #                 # #  #                  #     
                #               #                 # ###█###################     
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #                        
                #               #                 #    #     ################## 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               #                 #    #     #                # 
                #               ##################█#####     #                # 
                #                                 #          #                # 
                #                                 #          #                # 
                #                                 #          #                # 
                #                                 #          ################## 
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                #                                 #                             
                ###################################                             
""".strip('\n').split('\n')

print(find_rectangles(grid))