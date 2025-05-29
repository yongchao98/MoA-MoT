def find_rectangles(grid_str):
    # Split the grid into lines and remove any empty lines
    grid = [list(line) for line in grid_str.splitlines() if line]
    
    if not grid:
        return 0
        
    height = len(grid)
    width = len(grid[0])
    rectangles = 0

    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check horizontal lines (top and bottom)
        for x in range(left, right + 1):
            if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return False
        
        # Check vertical lines (left and right)
        for y in range(top, bottom + 1):
            if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return False
        return True

    # Find all potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] in ['#', '█']:
                # Found a potential top-left corner
                for y2 in range(y1 + 1, height):
                    if grid[y2][x1] in ['#', '█']:  # Found potential bottom
                        for x2 in range(x1 + 1, width):
                            if grid[y1][x2] in ['#', '█'] and grid[y2][x2] in ['#', '█']:
                                # Found potential rectangle, verify it
                                if is_rectangle(y1, x1, y2, x2):
                                    rectangles += 1
                                    break  # Move to next potential bottom edge

    return rectangles

# Create the test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                   #######################                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #                     #                      
                                   #######################                      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                          ######
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                 ##############################                           #    #
                 #                            #                           #    #
                 #                            #                           #    #
                 ##############################                           ######
                                                                                
                                                                                
                                                                                
                                #####################                           
                                #                   #                           
                                #####################                           
                                                                                """

print(find_rectangles(grid))