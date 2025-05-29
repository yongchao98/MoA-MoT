def find_rectangles(grid):
    # Convert input string to list of lines
    lines = grid.split('\n')
    rows = len(lines)
    cols = len(lines[0]) if rows > 0 else 0
    count = 0
    
    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check if all corners are present
        if (lines[top][left] not in '#█' or 
            lines[top][right] not in '#█' or 
            lines[bottom][left] not in '#█' or 
            lines[bottom][right] not in '#█'):
            return False
            
        # Check top and bottom edges
        for c in range(left + 1, right):
            if lines[top][c] not in '#█' or lines[bottom][c] not in '#█':
                return False
                
        # Check left and right edges
        for r in range(top + 1, bottom):
            if lines[r][left] not in '#█' or lines[r][right] not in '#█':
                return False
                
        # Check if interior is empty (or contains spaces)
        for r in range(top + 1, bottom):
            for c in range(left + 1, right):
                if lines[r][c] not in ' █':
                    return False
        
        return True
    
    # Search for potential rectangles
    for top in range(rows):
        for left in range(cols):
            # If we find a top-left corner
            if lines[top][left] in '#█':
                # Look for matching bottom-right corners
                for bottom in range(top + 1, rows):
                    for right in range(left + 1, cols):
                        if is_rectangle(top, left, bottom, right):
                            count += 1
    
    return count

# Test with the provided grid
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
                                    ############                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))