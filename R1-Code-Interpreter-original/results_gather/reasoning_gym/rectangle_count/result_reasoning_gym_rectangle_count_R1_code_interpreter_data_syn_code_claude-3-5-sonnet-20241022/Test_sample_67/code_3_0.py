def find_rectangles(grid):
    # Convert grid to list of strings
    lines = grid.split('\n')
    # Remove empty lines and ensure all lines have same width
    lines = [line for line in lines if line.strip()]
    
    def is_rectangle(top, left, bottom, right):
        # Check if the coordinates form a valid rectangle
        
        # Check top and bottom edges
        top_edge = lines[top][left:right+1]
        bottom_edge = lines[bottom][left:right+1]
        if not all(c in '#█' for c in top_edge) or not all(c in '#█' for c in bottom_edge):
            return False
            
        # Check left and right edges
        for row in range(top, bottom + 1):
            if lines[row][left] not in '#█' or lines[row][right] not in '#█':
                return False
                
        return True
    
    rectangles = set()
    height = len(lines)
    width = len(lines[0])
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            # Skip if not a potential corner
            if lines[top][left] not in '#█':
                continue
                
            # Look for bottom-right corners
            for bottom in range(top, height):
                for right in range(left, width):
                    if right > left and bottom > top:  # Ensure non-zero size
                        if is_rectangle(top, left, bottom, right):
                            # Verify it's a proper rectangle (not part of a larger one)
                            is_proper = False
                            # Check if it has at least one '#' character
                            for i in range(top, bottom + 1):
                                for j in range(left, right + 1):
                                    if lines[i][j] == '#':
                                        is_proper = True
                                        break
                                if is_proper:
                                    break
                            
                            if is_proper:
                                rectangles.add((top, left, bottom, right))
    
    # Count distinct rectangles
    count = 0
    processed = set()
    
    for rect in rectangles:
        top, left, bottom, right = rect
        
        # Skip if we've processed this area before
        if rect in processed:
            continue
            
        # Mark as processed
        processed.add(rect)
        
        # Check if it's a valid distinct rectangle
        is_distinct = True
        for other in rectangles:
            if other != rect and other not in processed:
                o_top, o_left, o_bottom, o_right = other
                # Check if completely overlapping
                if (o_top <= top and o_bottom >= bottom and 
                    o_left <= left and o_right >= right):
                    is_distinct = False
                    break
        
        if is_distinct:
            count += 1
    
    return count

# Test with the grid
test_grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                             ##################                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
                             #                #                                 
           ##################█##############################################    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #                 #               ###                           #    
           #               ##█#########      ###                           #    
           #               # #        #      ###                           #    
           #               # #        #      ###                           #    
           #               # #        #      ###                    ####   #    
           #               # #        #      ###                    #  #   #    
           #               # #        #      ###                    #  #   #    
           #               # #        #      ###                    #  #   #    
           #               ##█#########      ###                    #  #   #    
           ##################█████████████████##                    #  #   #    
                                               #                    #  #   #    
                                               #                    #  #   #    
                                               #                    ####   #    
                                               #                           #    
                                               #                           #    
                                               #                           #    
                                               #                           #    
                                               #     ######################█### 
                                               ######█###                  #  # 
                                               ##    ###█##################█### 
                                               #█████████###################    
                                                                                """

print(find_rectangles(test_grid))