def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    height = len(lines)
    width = len(lines[0]) if height > 0 else 0
    
    def find_rectangle_bounds(start_y, start_x):
        # Find right bound
        right = start_x
        while right < width and lines[start_y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom bound
        bottom = start_y
        while bottom < height and lines[bottom][start_x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        for y in range(start_y, bottom + 1):
            if lines[y][start_x] not in '#█' or lines[y][right] not in '#█':
                return None
        for x in range(start_x, right + 1):
            if lines[start_y][x] not in '#█' or lines[bottom][x] not in '#█':
                return None
                
        return (start_y, start_x, bottom, right)
    
    rectangles = set()
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if lines[y][x] in '#█':
                # Check if this is a top-left corner
                if (y == 0 or lines[y-1][x] not in '#█') and (x == 0 or lines[y][x-1] not in '#█'):
                    rect = find_rectangle_bounds(y, x)
                    if rect:
                        rectangles.add(rect)
    
    # The overlapping rectangles are already counted separately
    # because we're finding complete rectangles from their top-left corners
    print(len(rectangles))

# Create test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                 ############                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #             ###   
                                                 #          #             # #   
                                                 #          #             # #   
                                             ####█########  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             ####█████████###             # #   
                                                                          # #   
                                                                          # #   
                                                                          ###   
                                                                                
                                                                                """

find_rectangles(grid)