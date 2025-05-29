def find_rectangles(grid):
    lines = [line for line in grid.splitlines() if line.strip()]
    height = len(lines)
    width = len(lines[0]) if lines else 0
    
    def is_valid_rectangle(y1, x1, y2, x2):
        # Check if corners are valid
        if not (0 <= y1 < y2 < height and 0 <= x1 < x2 < width):
            return False
            
        # Check if it's a proper rectangle
        for y in range(y1, y2 + 1):
            for x in range(x1, x2 + 1):
                if y in (y1, y2) or x in (x1, x2):  # Only check the border
                    if lines[y][x] not in '#█':
                        return False
                elif lines[y][x] not in ' █':  # Inside should be empty or overlap
                    return False
        return True
    
    def has_overlap(y1, x1, y2, x2):
        return any(lines[y][x] == '█' 
                  for y in range(y1, y2 + 1)
                  for x in range(x1, x2 + 1))
    
    rectangles = set()
    # Find all potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if lines[y1][x1] in '#█':
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if lines[y2][x2] in '#█':
                            if is_valid_rectangle(y1, x1, y2, x2):
                                # Store rectangle coordinates
                                rectangles.add((y1, x1, y2, x2))
    
    # Count rectangles including overlaps
    total = 0
    for rect in rectangles:
        y1, x1, y2, x2 = rect
        # Add 2 if there's an overlap, 1 if not
        if has_overlap(y1, x1, y2, x2):
            total += 2
        else:
            total += 1
            
    return len(rectangles)  # Each distinct rectangle shape counts as 1

# Test with the provided ASCII art
ascii_art = """                                                                                
                                                                                
                                      ########################                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      ########################                  
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
                                                             # #                
                                                             # #                
                                                             # #                
                                                     ########█#█#######         
                                                     #       # #      #         
        ####################                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         ########█#█#######         
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        ####################                                 # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             ###                
                                                                                
                           ###################################################  
                           #                                                 #  
                           #                                                 #  
                           #                                                 #  
                           ###################################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(ascii_art))