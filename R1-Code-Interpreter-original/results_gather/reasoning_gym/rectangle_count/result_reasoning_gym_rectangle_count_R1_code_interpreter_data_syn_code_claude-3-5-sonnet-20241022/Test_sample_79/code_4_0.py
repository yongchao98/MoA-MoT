def find_rectangles(grid):
    lines = [line for line in grid.splitlines() if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def is_corner(x, y):
        if lines[y][x] not in '#█':
            return False
        
        # Count adjacent border cells
        adjacent = 0
        if y > 0 and lines[y-1][x] in '#█': adjacent += 1
        if y < height-1 and lines[y+1][x] in '#█': adjacent += 1
        if x > 0 and lines[y][x-1] in '#█': adjacent += 1
        if x < width-1 and lines[y][x+1] in '#█': adjacent += 1
        
        return adjacent == 2

    def verify_rectangle(x1, y1, x2, y2):
        if x1 >= x2 or y1 >= y2:
            return False
            
        # Check all corners
        if not all(lines[y][x] in '#█' for x, y in [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]):
            return False
            
        # Check horizontal borders
        for x in range(x1, x2 + 1):
            if lines[y1][x] not in '#█' or lines[y2][x] not in '#█':
                return False
                
        # Check vertical borders
        for y in range(y1, y2 + 1):
            if lines[y][x1] not in '#█' or lines[y][x2] not in '#█':
                return False
        
        return True

    # Find all corners
    corners = [(x, y) for y in range(height) for x in range(width) if is_corner(x, y)]
    
    rectangles = set()
    # Try to form rectangles from corners
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            # Try both possible rectangle configurations
            if verify_rectangle(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)):
                rectangles.add((min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)))
            
            # For overlapping rectangles, check alternative corners
            if lines[y1][x2] in '#█' and lines[y2][x1] in '#█':
                if verify_rectangle(x1, y1, x2, y2):
                    rectangles.add((x1, y1, x2, y2))

    # Additional check for overlapping rectangles
    for rect1 in rectangles.copy():
        x1, y1, x2, y2 = rect1
        # Check if there's an overlap point (█) that might form another rectangle
        for y in range(y1, y2 + 1):
            for x in range(x1, x2 + 1):
                if lines[y][x] == '█':
                    # Check if this point forms part of another rectangle
                    for rect2 in rectangles.copy():
                        if rect2 != rect1:
                            rx1, ry1, rx2, ry2 = rect2
                            if rx1 <= x <= rx2 and ry1 <= y <= ry2:
                                # Check if overlap forms a new rectangle
                                new_corners = [(rx1, y1), (rx2, y2), (x1, ry1), (x2, ry2)]
                                for c1 in range(4):
                                    for c2 in range(c1 + 1, 4):
                                        nx1, ny1 = new_corners[c1]
                                        nx2, ny2 = new_corners[c2]
                                        if verify_rectangle(min(nx1, nx2), min(ny1, ny2), 
                                                          max(nx1, nx2), max(ny1, ny2)):
                                            rectangles.add((min(nx1, nx2), min(ny1, ny2), 
                                                          max(nx1, nx2), max(ny1, ny2)))

    return len(rectangles)

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