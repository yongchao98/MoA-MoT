def is_valid_corner(grid, x, y):
    return 0 <= y < len(grid) and 0 <= x < len(grid[0]) and grid[y][x] in '#█'

def verify_rectangle(grid, x1, y1, x2, y2):
    if x1 >= x2 or y1 >= y2:
        return False
    
    # Check all four corners
    if not all(is_valid_corner(grid, x, y) for x, y in [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]):
        return False
    
    # Check top and bottom edges
    for x in range(x1 + 1, x2):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
    
    # Check left and right edges
    for y in range(y1 + 1, y2):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
    
    return True

def find_all_rectangles(grid):
    # Clean and prepare the grid
    lines = [line for line in grid.split('\n') if line.strip()]
    grid = [list(line) for line in lines]
    height = len(grid)
    width = len(grid[0])
    
    # Find all corner points
    corner_points = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                corner_points.append((x, y))
    
    rectangles = set()
    # Check all possible combinations of corner points
    for i, (x1, y1) in enumerate(corner_points):
        for j, (x2, y2) in enumerate(corner_points[i+1:], i+1):
            # Try to form rectangles with these points as opposite corners
            if x1 < x2 and y1 < y2:
                if verify_rectangle(grid, x1, y1, x2, y2):
                    rectangles.add((x1, y1, x2, y2))
            
            # Also try other corner combinations
            if x2 < x1 and y1 < y2:
                if verify_rectangle(grid, x2, y1, x1, y2):
                    rectangles.add((x2, y1, x1, y2))
    
    return len(rectangles)

# Input grid
grid = """   ###############                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #             #                                                              
   #          ###█#####################################################         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          #  #                                                    #         
   #          ###█##█████████##############################█████#######         
   #             #  #       █##################            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #   #                
   #             #  #       █                 #            #####                
   #             #  #       █                 #                                 
   #             #  #       █                 #                                 
   #             #  #       █                 #                                 
   ###############  #       █                 #                                 
                    #       █                 #                                 
                    ########█                 #                                 
                            #                 #                                 
                            #                 #                                 
                           #█#############    #                                 
                           ##            #    #                                 
                           ##############█#####                                 
                           #             #                                      
                           #             #                                      
                           #             #       #############################  
                           #             #       #                           #  
                           #             #       #                           #  
                           #             #       #                           #  
                           #             #       #                           #  
                           ###############       #                           #  
                               ###############   #                           #  
                               #             #   #                           #  
                               #             #   ███████████████##############  
                               #             #   #             #                
                               ###############   ###############                """

print(find_all_rectangles(grid))