def find_corners(grid):
    corners = []
    height = len(grid)
    width = len(grid[0])
    
    # Check each point in the grid
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Count adjacent '#' or '█' characters
                neighbors = 0
                if y > 0 and grid[y-1][x] in '#█': neighbors += 1
                if y < height-1 and grid[y+1][x] in '#█': neighbors += 1
                if x > 0 and grid[y][x-1] in '#█': neighbors += 1
                if x < width-1 and grid[y][x+1] in '#█': neighbors += 1
                
                # If point has exactly 2 adjacent '#' or '█', it's a corner
                if neighbors == 2:
                    corners.append((x, y))
    return corners

def is_valid_rectangle(corners, grid):
    if len(corners) != 4:
        return False
    
    # Sort corners by x, then y
    corners = sorted(corners)
    
    # Check if corners form a rectangle
    x1, y1 = corners[0]
    x2, y2 = corners[1]
    x3, y3 = corners[2]
    x4, y4 = corners[3]
    
    # Verify edges
    for y in range(min(y1, y2), max(y1, y2) + 1):
        if grid[y][x1] not in '#█':
            return False
    for y in range(min(y3, y4), max(y3, y4) + 1):
        if grid[y][x3] not in '#█':
            return False
    for x in range(min(x1, x3), max(x1, x3) + 1):
        if grid[y1][x] not in '#█':
            return False
    for x in range(min(x2, x4), max(x2, x4) + 1):
        if grid[y4][x] not in '#█':
            return False
    
    return True

def count_rectangles(grid):
    # Convert input string to list of lines
    grid = grid.splitlines()
    
    # Find all corners
    corners = find_corners(grid)
    
    # Try all possible combinations of 4 corners
    rectangle_count = 0
    from itertools import combinations
    
    for corner_set in combinations(corners, 4):
        if is_valid_rectangle(corner_set, grid):
            rectangle_count += 1
    
    return rectangle_count

# Create the grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                     ######################     
                                                     #                    #     
                                                     #                    #     
                                        #############█#####################     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #############█#####################     
                                                     #                    #     
                        #############################█#################   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #    ##############          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #################█####     
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #####█############█############################         
                             #            #                                     
                             ##############                                     
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(count_rectangles(grid))