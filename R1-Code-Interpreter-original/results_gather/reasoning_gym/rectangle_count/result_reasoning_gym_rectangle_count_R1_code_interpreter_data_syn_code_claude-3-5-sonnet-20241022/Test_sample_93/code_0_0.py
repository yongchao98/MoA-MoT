def is_corner(grid, i, j):
    # Check if position is a corner (has perpendicular edges)
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    if grid[i][j] not in '#█':
        return False
    
    # Count adjacent walls
    adjacent = 0
    horizontal = 0
    vertical = 0
    
    # Check horizontal
    if j > 0 and grid[i][j-1] in '#█':
        horizontal += 1
    if j < len(grid[0])-1 and grid[i][j+1] in '#█':
        horizontal += 1
        
    # Check vertical
    if i > 0 and grid[i-1][j] in '#█':
        vertical += 1
    if i < len(grid)-1 and grid[i+1][j] in '#█':
        vertical += 1
        
    return (horizontal == 1 and vertical == 1)

def find_rectangle(grid, start_i, start_j):
    # Find complete rectangle starting from top-left corner
    if not is_corner(grid, start_i, start_j):
        return None
    
    # Find right edge
    right_j = start_j
    while right_j < len(grid[0]) and grid[start_i][right_j] in '#█':
        right_j += 1
    right_j -= 1
    
    # Find bottom edge
    bottom_i = start_i
    while bottom_i < len(grid) and grid[bottom_i][start_j] in '#█':
        bottom_i += 1
    bottom_i -= 1
    
    # Validate rectangle
    if not (is_corner(grid, start_i, right_j) and 
            is_corner(grid, bottom_i, start_j) and 
            is_corner(grid, bottom_i, right_j)):
        return None
    
    # Check if all edges are complete
    for i in range(start_i, bottom_i+1):
        if grid[i][start_j] not in '#█' or grid[i][right_j] not in '#█':
            return None
    for j in range(start_j, right_j+1):
        if grid[start_i][j] not in '#█' or grid[bottom_i][j] not in '#█':
            return None
            
    return (start_i, start_j, bottom_i, right_j)

def count_rectangles(grid):
    rectangles = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in '#█':
                rect = find_rectangle(grid, i, j)
                if rect:
                    rectangles.add(rect)
    return len(rectangles)

# Create grid from the ASCII art
grid = []
empty_lines = 0
for line in """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                       ###                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                                       # #                                      
                       ################█#█######                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               # #     #                                
                       #               ###     #                                
                       #                       #                                
                       #                       #                                
                       #                       #                                
                       #                       #                                
                       #                       #                                
                       #                       #                                
                       #                       #                                
                       #    ###################█#####################           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #                    #           
                       #    #                  #            ########█########## 
                       #    #                  #            #       #         # 
                       #    #                  #            #       #         # 
                       #    #                  #            #       #         # 
                       #    #                  #            #       #         # 
                       #    #                  #            #       #         # 
                       #####█###################            #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               #       #         # 
                            #                               ########█########## 
                            #########################################           
                                                                                
                                                                                
                                                                                
""".splitlines():
    if line.strip():
        grid.append(list(line))
    
print(count_rectangles(grid))