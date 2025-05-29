def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def find_rectangle_end(start_x, start_y):
        # Find width of top line
        width = 0
        for x in range(start_x, width):
            if grid[start_y][x] != '#':
                break
            width += 1
        if width < 2:  # Need at least 2 wide
            return None
            
        # Find height by checking right edge
        height = 0
        for y in range(start_y, height):
            if grid[y][start_x + width - 1] != '#':
                break
            height += 1
        if height < 2:  # Need at least 2 high
            return None
            
        # Verify left edge
        for y in range(start_y, start_y + height):
            if grid[y][start_x] != '#':
                return None
                
        # Verify bottom edge
        for x in range(start_x, start_x + width):
            if grid[start_y + height - 1][x] != '#':
                return None
                
        return (width, height)
    
    rectangles = 0
    # Scan for top-left corners (must have right and bottom segments)
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '#':
                # Check if it's a top-left corner
                if (x == 0 or grid[y][x-1] != '#') and (y == 0 or grid[y-1][x] != '#'):
                    if x + 1 < width and y + 1 < height:
                        if grid[y][x+1] == '#' and grid[y+1][x] == '#':
                            rectangles += 1
    
    print(rectangles)

# Create the grid string (with proper handling of the last empty line)
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                      #######                                   
                                      #     #                                   
                                      #     #                                   
                                      #     #                                   
                                      #     #                                   
                                      #######                ###                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             ###                
                                                                                
                                                                                
                                                                                
                                                                                
                                                   ###################          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   ###################          
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                     #######    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #######"""

find_rectangles(grid)