def find_rectangles(grid):
    # Convert input to lines and filter out empty lines
    lines = [line for line in grid.splitlines() if line.strip()]
    
    if not lines:
        return 0
    
    # Find all '#' positions
    points = []
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if lines[y][x] == '#':
                points.append((y, x))
    
    def is_minimal_rectangle(y1, x1, y2, x2):
        # Check if there are any '#' characters inside the rectangle
        # (excluding the borders)
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if lines[y][x] == '#':
                    # If there's a '#' that's not part of the edges,
                    # this is not a minimal rectangle
                    if not (y in [y1, y2] or x in [x1, x2]):
                        return False
        return True
    
    # Find rectangles by checking all possible corners
    rectangles = set()
    for i, (y1, x1) in enumerate(points):
        for y2, x2 in points[i:]:
            # Ensure points form a diagonal of a rectangle (y2 > y1 and x2 > x1)
            if y2 <= y1 or x2 <= x1:
                continue
                
            # Check if other two corners exist
            if (y1, x2) not in points or (y2, x1) not in points:
                continue
                
            # Verify edges
            valid = True
            # Check horizontal edges
            for x in range(x1 + 1, x2):
                if lines[y1][x] != '#' or lines[y2][x] != '#':
                    valid = False
                    break
                    
            if not valid:
                continue
                
            # Check vertical edges
            for y in range(y1 + 1, y2):
                if lines[y][x1] != '#' or lines[y][x2] != '#':
                    valid = False
                    break
                    
            # Only add if it's a valid rectangle and is minimal
            if valid and is_minimal_rectangle(y1, x1, y2, x2):
                rectangles.add((y1, x1, y2, x2))
    
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