def is_complete_rectangle(grid, top, left, bottom, right):
    # Check if there are any breaks in the rectangle outline
    # and ensure the inside is empty (except for overlaps)
    
    # Check horizontal edges
    for x in range(left, right + 1):
        # Top edge must be continuous
        if grid[top][x] not in '#█':
            return False
        # Bottom edge must be continuous
        if grid[bottom][x] not in '#█':
            return False
        
    # Check vertical edges
    for y in range(top, bottom + 1):
        # Left edge must be continuous
        if grid[y][left] not in '#█':
            return False
        # Right edge must be continuous
        if grid[y][right] not in '#█':
            return False
        
        # Check inside of rectangle
        for x in range(left + 1, right):
            if y != top and y != bottom:  # Skip edges
                if grid[y][x] not in ' #█':  # Inside should be empty or part of another rectangle
                    return False
    
    # Check minimum size
    if right - left < 2 or bottom - top < 2:
        return False
        
    return True

def find_rectangles(grid):
    # Split the grid into lines and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    # Convert to 2D array
    grid = [list(row) for row in lines]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = set()
    
    # Find all possible rectangles
    for top in range(height - 2):  # Need at least 2 rows
        for left in range(width - 2):  # Need at least 2 columns
            # Only start from corner pieces
            if grid[top][left] not in '#█':
                continue
            
            for bottom in range(top + 2, height):  # Minimum height of 2
                for right in range(left + 2, width):  # Minimum width of 2
                    if grid[bottom][right] not in '#█':
                        continue
                        
                    if is_complete_rectangle(grid, top, left, bottom, right):
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

print(find_rectangles(grid))