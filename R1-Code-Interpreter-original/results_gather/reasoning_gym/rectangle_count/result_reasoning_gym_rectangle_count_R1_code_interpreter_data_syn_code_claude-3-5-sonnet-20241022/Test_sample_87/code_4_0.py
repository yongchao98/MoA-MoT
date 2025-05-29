def find_rectangles(grid_str):
    # Convert grid string to list of lines and remove any trailing spaces
    grid = [line.rstrip() for line in grid_str.strip().split('\n')]
    height = len(grid)
    width = len(grid[0])
    
    def find_horizontal_line(y, x1, x2):
        return all(grid[y][x] in '#█' for x in range(x1, x2 + 1))
    
    def find_vertical_line(x, y1, y2):
        return all(grid[y][x] in '#█' for y in range(y1, y2 + 1))
    
    rectangles = set()
    
    # Find all potential top-left corners
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] not in '#█':
                continue
                
            # Look for right edge of rectangle
            for x2 in range(x1 + 1, width):
                if grid[y1][x2] not in '#█':
                    continue
                    
                # Check if we have a horizontal line at the top
                if not find_horizontal_line(y1, x1, x2):
                    continue
                
                # Look for bottom edge
                for y2 in range(y1 + 1, height):
                    if grid[y2][x1] not in '#█' or grid[y2][x2] not in '#█':
                        continue
                        
                    # Check if we have a complete rectangle
                    if (find_horizontal_line(y2, x1, x2) and
                        find_vertical_line(x1, y1, y2) and
                        find_vertical_line(x2, y1, y2)):
                        
                        # Check if interior is valid (either empty or contains █)
                        valid = True
                        for y in range(y1 + 1, y2):
                            for x in range(x1 + 1, x2):
                                if grid[y][x] == '#':
                                    valid = False
                                    break
                            if not valid:
                                break
                                
                        if valid:
                            rectangles.add((y1, x1, y2, x2))
    
    return len(rectangles)

# Define the test grid (same as before)
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                       ######                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                               ##############      ####█####█###############    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               ##############      #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   ####█####█###############    
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       ######                   
                                                                                
                       ##################################################       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #  #################################             #       
                       #  #                               #             #       
                       #  #################################             #       
                       #                                                #       
                       ##################################################       
                                                                                """

print(find_rectangles(grid))