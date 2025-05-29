def find_rectangles(grid):
    # Convert string grid to list of lines
    lines = grid.split('\n')
    
    # Find all '#' positions
    hash_positions = []
    for y in range(len(lines)):
        for x in range(len(lines[y])):
            if lines[y][x] in ['#', '█']:
                hash_positions.append((x, y))
    
    # Find rectangles
    rectangles = set()
    for x1, y1 in hash_positions:
        # Look for potential top-left corners
        if (x1-1, y1) not in hash_positions and (x1, y1-1) not in hash_positions:
            # For each potential top-left corner, look for matching corners
            for x2, y2 in hash_positions:
                if x2 > x1 and y2 > y1:  # potential bottom-right corner
                    # Check if it forms a rectangle
                    is_rectangle = True
                    # Check top and bottom edges
                    for x in range(x1, x2+1):
                        if (x, y1) not in hash_positions or (x, y2) not in hash_positions:
                            is_rectangle = False
                            break
                    # Check left and right edges
                    for y in range(y1, y2+1):
                        if (x1, y) not in hash_positions or (x2, y) not in hash_positions:
                            is_rectangle = False
                            break
                    if is_rectangle:
                        rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Create the grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                       ##############################           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       #                            #           
                                       ##############################           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                     ###################                        
                                     #                 #                        
                                     ###################                        
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))