def count_rectangles(grid):
    # Convert the grid into a list of strings
    lines = grid.splitlines()
    # Remove any empty lines
    lines = [line for line in lines if line]
    
    if not lines:
        return 0
        
    height = len(lines)
    width = len(lines[0])
    
    # Find all '#' positions
    corners = []
    for y in range(height):
        for x in range(width):
            if lines[y][x] == '#':
                corners.append((x, y))
    
    if not corners:
        return 0
    
    # Find rectangles
    rectangles = 0
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            # Check if these points could form a rectangle
            if x1 < x2 and y1 < y2:
                # Check if other two corners exist
                if (x1, y2) in corners and (x2, y1) in corners:
                    # Verify all connecting lines exist
                    is_rectangle = True
                    # Check horizontal lines
                    for x in range(x1+1, x2):
                        if (x, y1) not in corners or (x, y2) not in corners:
                            is_rectangle = False
                            break
                    if not is_rectangle:
                        continue
                    # Check vertical lines
                    for y in range(y1+1, y2):
                        if (x1, y) not in corners or (x2, y) not in corners:
                            is_rectangle = False
                            break
                    if is_rectangle:
                        rectangles += 1
    
    return rectangles

# Create test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                       #########
                                                                       #       #
                                                                       #########
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(count_rectangles(grid))