def find_rectangles(grid):
    # Convert grid to list of lines and remove empty lines
    lines = [line for line in grid.splitlines() if line.strip()]
    
    # Find all vertical lines (sequences of # or █)
    verticals = []
    for x in range(len(lines[0])):
        start = None
        for y in range(len(lines)):
            if lines[y][x] in '#█':
                if start is None:
                    start = y
            elif start is not None:
                verticals.append((x, start, y-1))
                start = None
        if start is not None:
            verticals.append((x, start, len(lines)-1))

    # Find all horizontal lines
    horizontals = []
    for y in range(len(lines)):
        start = None
        for x in range(len(lines[0])):
            if lines[y][x] in '#█':
                if start is None:
                    start = x
            elif start is not None:
                horizontals.append((y, start, x-1))
                start = None
        if start is not None:
            horizontals.append((y, start, len(lines[0])-1))

    # Count rectangles by matching horizontal and vertical lines
    rectangles = 0
    for v1 in verticals:
        for v2 in verticals:
            if v1[0] >= v2[0]:
                continue
            # For each pair of vertical lines
            x1, y1_start, y1_end = v1
            x2, y2_start, y2_end = v2
            
            # Look for matching horizontal lines
            for h1 in horizontals:
                if h1[1] > x1 or h1[2] < x2:
                    continue
                y, x_start, x_end = h1
                
                # Check for second horizontal line
                for h2 in horizontals:
                    if h2[1] > x1 or h2[2] < x2 or h2[0] <= y:
                        continue
                    y2, x2_start, x2_end = h2
                    
                    # Verify rectangle
                    if (x_start <= x1 and x_end >= x2 and 
                        x2_start <= x1 and x2_end >= x2 and
                        y >= y1_start and y <= y1_end and
                        y2 >= y2_start and y2 <= y2_end):
                        rectangles += 1

    return rectangles

# Create the grid (same as before)
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                     ######################     
                                                     #                    #     
                                                     #                    #     
                                        #############█#####################     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #            #                   ##     
                                        #############█#####################     
                                                     #                    #     
                        #############################█#################   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #                            #                #   #     
                        #    ##############          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #                #   #     
                        #    #            #          #################█####     
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #    #            #                           #         
                        #####█############█############################         
                             #            #                                     
                             ##############                                     
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))