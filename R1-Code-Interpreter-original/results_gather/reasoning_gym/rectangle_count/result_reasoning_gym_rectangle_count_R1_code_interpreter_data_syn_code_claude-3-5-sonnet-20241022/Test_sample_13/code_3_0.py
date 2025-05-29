def find_rectangles(grid):
    lines = [line for line in grid.splitlines() if line.strip()]
    
    def is_complete_rectangle(y, x):
        if lines[y][x] not in '#█':
            return False
            
        # Find right edge
        right = x + 1
        while right < len(lines[y]) and lines[y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = y + 1
        while bottom < len(lines) and lines[bottom][x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify it's a complete rectangle
        if lines[y][right] not in '#█' or lines[bottom][x] not in '#█' or lines[bottom][right] not in '#█':
            return False
            
        # Check all borders are complete
        for i in range(x, right + 1):
            if lines[y][i] not in '#█' or lines[bottom][i] not in '#█':
                return False
        for i in range(y, bottom + 1):
            if lines[i][x] not in '#█' or lines[i][right] not in '#█':
                return False
                
        # Check if there's an overlap
        has_overlap = False
        for i in range(y, bottom + 1):
            for j in range(x, right + 1):
                if lines[i][j] == '█':
                    has_overlap = True
                    
        return 2 if has_overlap else 1
    
    # Count distinct rectangles
    count = 0
    counted_positions = set()
    
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if (y, x) not in counted_positions and lines[y][x] in '#█':
                result = is_complete_rectangle(y, x)
                if result:
                    count += result
                    # Mark this position as counted
                    counted_positions.add((y, x))
    
    return count

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