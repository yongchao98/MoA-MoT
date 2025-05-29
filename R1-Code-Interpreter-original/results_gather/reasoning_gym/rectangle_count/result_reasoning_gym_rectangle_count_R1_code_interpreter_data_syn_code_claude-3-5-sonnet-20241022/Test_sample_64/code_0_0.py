def is_valid_rectangle(grid, top, left, bottom, right):
    # Check if corners are marked
    corners = [grid[top][left], grid[top][right], grid[bottom][left], grid[bottom][right]]
    if not all(c in '#█' for c in corners):
        return False
    
    # Check horizontal edges
    for x in range(left + 1, right):
        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
            return False
    
    # Check vertical edges
    for y in range(top + 1, bottom):
        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
            return False
    
    return True

def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid]
    
    # Find all rectangles
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in '#█':
                continue
            
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if is_valid_rectangle(grid, top, left, bottom, right):
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

print(count_rectangles(grid))