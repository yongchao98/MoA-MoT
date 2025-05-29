def find_rectangles(grid):
    # Convert string to list of lines and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    height = len(lines)
    width = len(lines[0]) if lines else 0
    rectangles = 0
    
    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        if not all(lines[top][x] in '#█' for x in range(left, right + 1)):
            return False
        if not all(lines[bottom][x] in '#█' for x in range(left, right + 1)):
            return False
        
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not lines[y][left] in '#█' or not lines[y][right] in '#█':
                return False
        return True
    
    # Scan for potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if lines[y1][x1] in '#█':
                # Found potential top-left corner
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if lines[y2][x2] in '#█':
                            # Check if this forms a valid rectangle
                            if is_rectangle(y1, x1, y2, x2):
                                rectangles += 1
                                # Skip to next potential start point to avoid counting the same rectangle
                                break
                    
    return rectangles

# Create the grid string (your provided grid)
grid = """
                                                                                
                                                                                
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
                                                               ###              
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
                           #          #                                         
                           #          #                                         
                           #          #                                         
                           ############                                         
"""

print(find_rectangles(grid))