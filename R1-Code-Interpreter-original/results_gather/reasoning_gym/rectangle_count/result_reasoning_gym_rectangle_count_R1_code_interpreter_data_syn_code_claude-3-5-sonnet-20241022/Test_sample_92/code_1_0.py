def find_rectangles(grid):
    # Convert input to lines
    lines = grid.splitlines()
    
    # Find all continuous horizontal and vertical lines
    horizontals = []  # (y, x_start, x_end)
    verticals = []    # (x, y_start, y_end)
    
    # Find horizontal lines
    for y in range(len(lines)):
        start = None
        for x in range(len(lines[y])):
            if lines[y][x] == '#':
                if start is None:
                    start = x
            elif start is not None:
                if x - start > 1:  # Minimum length of 2 to be part of rectangle
                    horizontals.append((y, start, x-1))
                start = None
        if start is not None and len(lines[y]) - start > 1:
            horizontals.append((y, start, len(lines[y])-1))
    
    # Find vertical lines
    for x in range(len(lines[0])):
        start = None
        for y in range(len(lines)):
            if lines[y][x] == '#':
                if start is None:
                    start = y
            elif start is not None:
                if y - start > 1:  # Minimum length of 2 to be part of rectangle
                    verticals.append((x, start, y-1))
                start = None
        if start is not None and len(lines) - start > 1:
            verticals.append((x, start, len(lines)-1))
    
    # Find rectangles by matching horizontal and vertical lines
    rectangles = set()
    for h1 in horizontals:
        for h2 in horizontals:
            if h1[0] >= h2[0]:  # Ensure h1 is above h2
                continue
            # Check if the horizontal lines have the same width
            if h1[1] != h2[1] or h1[2] != h2[2]:
                continue
            # Look for matching vertical lines
            found_left = False
            found_right = False
            for v in verticals:
                if v[0] == h1[1]:  # Left vertical line
                    if v[1] == h1[0] and v[2] == h2[0]:
                        found_left = True
                elif v[0] == h1[2]:  # Right vertical line
                    if v[1] == h1[0] and v[2] == h2[0]:
                        found_right = True
            if found_left and found_right:
                rectangles.add((h1[0], h1[1], h2[0], h1[2]))
    
    return len(rectangles)

# Test grid
test_grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                      ###       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      # #       
                                                                      ###       
                                                                                
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
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      # #                       
                                                      ###                       
                                                                                
                                                     ###########                
                                                     #         #                
                                                     #         #                
                                                     #         #                
                                        ########################                
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #           #                           
                                        #############                           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                     ###                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     # #                        
                                                     ###                        
                                                                                """

print(find_rectangles(test_grid))