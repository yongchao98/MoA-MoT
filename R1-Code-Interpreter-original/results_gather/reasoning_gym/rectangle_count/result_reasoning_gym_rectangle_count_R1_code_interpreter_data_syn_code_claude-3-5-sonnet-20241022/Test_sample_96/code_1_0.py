def find_rectangles(grid):
    # Split the grid into lines and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    def is_top_left_corner(x, y):
        # Check if current position is '#'
        if lines[y][x] != '#':
            return False
        # Check if it's a top corner (no '#' above)
        if y > 0 and lines[y-1][x] == '#':
            return False
        # Check if it's a left corner (no '#' to the left)
        if x > 0 and lines[y][x-1] == '#':
            return False
        return True
    
    def get_rectangle_size(start_x, start_y):
        # Get width (follow top edge until no more '#')
        width = 0
        for x in range(start_x, len(lines[start_y])):
            if lines[start_y][x] != '#':
                break
            width += 1
        
        # Get height (follow left edge until no more '#')
        height = 0
        for y in range(start_y, len(lines)):
            if lines[y][start_x] != '#':
                break
            height += 1
        
        return width, height
    
    def verify_rectangle(x, y, width, height):
        if width < 2 or height < 2:
            return False
        
        # Check if we have enough space
        if y + height > len(lines) or x + width > len(lines[y]):
            return False
        
        # Verify right edge
        for dy in range(height):
            if lines[y + dy][x + width - 1] != '#':
                return False
        
        # Verify bottom edge
        for dx in range(width):
            if lines[y + height - 1][x + dx] != '#':
                return False
        
        return True
    
    count = 0
    # Scan through the grid
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if is_top_left_corner(x, y):
                width, height = get_rectangle_size(x, y)
                if verify_rectangle(x, y, width, height):
                    count += 1
    
    print(count)

# Create the grid string
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