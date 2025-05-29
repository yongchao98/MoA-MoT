def count_rectangles(grid):
    # Convert string to list of lines
    lines = grid.split('\n')
    
    # Find dimensions
    height = len(lines)
    width = len(lines[0]) if lines else 0
    
    # Initialize counter
    rectangle_count = 0
    
    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check horizontal lines
        for x in range(left, right + 1):
            if lines[top][x] not in '#█' or lines[bottom][x] not in '#█':
                return False
        # Check vertical lines
        for y in range(top, bottom + 1):
            if lines[y][left] not in '#█' or lines[y][right] not in '#█':
                return False
        return True
    
    # Scan for potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if lines[y1][x1] in '#█':
                # Found top-left corner, now look for bottom-right
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if lines[y2][x2] in '#█':
                            # Verify if it forms a complete rectangle
                            if is_rectangle(y1, x1, y2, x2):
                                # Verify it's not just a line
                                if y2 > y1 and x2 > x1:
                                    rectangle_count += 1
                                    # Skip to next potential top-left corner
                                    break
                    else:
                        continue
                    break
    
    print(rectangle_count)

# Your ASCII art grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                 ##################################################             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 #                                                #             
                 ##################################################             
                                                                                
                                                                                
                                                                                
                                                                                
   ######################################                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   #                                    #                                       
   ######################################                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

count_rectangles(grid)