def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    # Find all '#' and '█' positions
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Function to check if a point is part of a rectangle
    def is_valid_char(char):
        return char in '#█'
    
    # Function to trace rectangle from top-left corner
    def trace_rectangle(start_y, start_x):
        # Find width
        x = start_x
        while x < width and is_valid_char(grid[start_y][x]):
            x += 1
        rect_width = x - start_x
        
        # Find height
        y = start_y
        while y < height and is_valid_char(grid[y][start_x]):
            y += 1
        rect_height = y - start_y
        
        # Verify rectangle
        for cy in range(start_y, start_y + rect_height):
            for cx in range(start_x, start_x + rect_width):
                if not is_valid_char(grid[cy][cx]):
                    return None
        
        return (start_y, start_x, rect_height, rect_width)
    
    # Find all rectangles
    visited = set()
    for y in range(height):
        for x in range(width):
            if (y, x) not in visited and is_valid_char(grid[y][x]):
                rect = trace_rectangle(y, x)
                if rect:
                    rectangles.append(rect)
                    # Mark all points in rectangle as visited
                    for cy in range(rect[0], rect[0] + rect[2]):
                        for cx in range(rect[1], rect[1] + rect[3]):
                            visited.add((cy, cx))
    
    # Count overlapping areas
    overlap_count = sum(1 for y in range(height) 
                       for x in range(width) 
                       if grid[y][x] == '█')
    
    # Total rectangles = number of distinct rectangles + 2 * number of overlapping areas
    total = len(rectangles) + overlap_count
    
    print(total)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                    ####███#    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    #   # ##    
                                                                    ####█#█#    
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        # #     
                                                                        ###     
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
      ######################                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                                    
      #                    #                                      #######       
      #                    #                                      #     #       
      ######################                                      #######       
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)