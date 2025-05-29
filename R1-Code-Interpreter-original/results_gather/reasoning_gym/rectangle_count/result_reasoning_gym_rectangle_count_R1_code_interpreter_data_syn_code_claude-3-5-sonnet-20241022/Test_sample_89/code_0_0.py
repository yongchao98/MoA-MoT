def find_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [line for line in grid.split('\n') if line]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    # Find all '#' and '█' positions
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if this is a corner by looking at neighbors
                is_corner = False
                
                # Top-left corner
                if (x == 0 or grid[y][x-1] not in '#█') and (y == 0 or grid[y-1][x] not in '#█'):
                    if x < width-1 and y < height-1 and grid[y][x+1] in '#█' and grid[y+1][x] in '#█':
                        corners.append((x, y, 'TL'))
                
                # Top-right corner
                if (x == width-1 or grid[y][x+1] not in '#█') and (y == 0 or grid[y-1][x] not in '#█'):
                    if x > 0 and y < height-1 and grid[y][x-1] in '#█' and grid[y+1][x] in '#█':
                        corners.append((x, y, 'TR'))
                
                # Bottom-left corner
                if (x == 0 or grid[y][x-1] not in '#█') and (y == height-1 or grid[y+1][x] not in '#█'):
                    if x < width-1 and y > 0 and grid[y][x+1] in '#█' and grid[y-1][x] in '#█':
                        corners.append((x, y, 'BL'))
                
                # Bottom-right corner
                if (x == width-1 or grid[y][x+1] not in '#█') and (y == height-1 or grid[y+1][x] not in '#█'):
                    if x > 0 and y > 0 and grid[y][x-1] in '#█' and grid[y-1][x] in '#█':
                        corners.append((x, y, 'BR'))
    
    # Group corners into rectangles
    rectangles = []
    for corner in corners:
        x, y, corner_type = corner
        if corner_type == 'TL':
            # Find matching bottom-right corner
            for other in corners:
                ox, oy, other_type = other
                if other_type == 'BR' and ox > x and oy > y:
                    # Verify this forms a valid rectangle
                    is_valid = True
                    # Check top and bottom edges
                    for i in range(x, ox + 1):
                        if grid[y][i] not in '#█' or grid[oy][i] not in '#█':
                            is_valid = False
                            break
                    # Check left and right edges
                    for i in range(y, oy + 1):
                        if grid[i][x] not in '#█' or grid[i][ox] not in '#█':
                            is_valid = False
                            break
                    if is_valid:
                        rectangles.append((x, y, ox, oy))
    
    print(len(rectangles))

# Create the grid
grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                           ########################################             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                           #                                      #             
                      #####█####################################  #             
                      #    #                                   #  #             
                      #    #                                   #  #             
                      #    #                                   #  #             
                      #####█####################################  #             
                           #                                      #             
                           ########################################             
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

find_rectangles(grid)