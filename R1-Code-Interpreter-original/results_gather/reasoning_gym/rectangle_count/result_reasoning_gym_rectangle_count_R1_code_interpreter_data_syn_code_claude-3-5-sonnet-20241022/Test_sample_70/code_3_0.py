def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    height = len(lines)
    width = len(lines[0]) if height > 0 else 0
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if not (lines[top][left] in '#█' and lines[top][right] in '#█' and 
                lines[bottom][left] in '#█' and lines[bottom][right] in '#█'):
            return False
            
        # Check top and bottom edges
        for x in range(left, right + 1):
            if lines[top][x] not in '#█' or lines[bottom][x] not in '#█':
                return False
                
        # Check left and right edges
        for y in range(top, bottom + 1):
            if lines[y][left] not in '#█' or lines[y][right] not in '#█':
                return False
                
        # Check if there's content inside (must be either space or #)
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if lines[y][x] not in ' #█':
                    return False
        
        return True

    # Store unique rectangles
    rectangles = set()
    
    # Find the main rectangle (tall one in the middle)
    for i in range(height):
        for j in range(width):
            if lines[i][j] == '#':
                # Try to find complete rectangles starting from this point
                for b in range(i + 1, height):
                    for r in range(j + 1, width):
                        if is_valid_rectangle(i, j, b, r):
                            # Store rectangle coordinates
                            rectangles.add((i, j, b, r))
    
    # Count total rectangles
    total = len(rectangles)
    print(total)

# Create test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                 ############                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #                   
                                                 #          #             ###   
                                                 #          #             # #   
                                                 #          #             # #   
                                             ####█########  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             #   #       #  #             # #   
                                             ####█████████###             # #   
                                                                          # #   
                                                                          # #   
                                                                          ###   
                                                                                
                                                                                """

find_rectangles(grid)