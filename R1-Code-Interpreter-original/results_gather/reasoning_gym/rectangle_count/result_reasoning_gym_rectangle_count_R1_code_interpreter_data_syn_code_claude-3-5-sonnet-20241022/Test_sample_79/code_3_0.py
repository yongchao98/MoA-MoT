def find_rectangles(grid):
    lines = [line for line in grid.splitlines() if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def trace_rectangle(x, y):
        # Find width of top edge
        width = 0
        while x + width < len(lines[0]) and lines[y][x + width] in '#█':
            width += 1
        if width < 2:  # Too narrow to be a rectangle
            return None
            
        # Find height of left edge
        height = 0
        while y + height < len(lines) and lines[y + height][x] in '#█':
            height += 1
        if height < 2:  # Too short to be a rectangle
            return None
            
        # Verify right edge
        for h in range(height):
            if lines[y + h][x + width - 1] not in '#█':
                return None
                
        # Verify bottom edge
        for w in range(width):
            if lines[y + height - 1][x + w] not in '#█':
                return None
                
        # Verify interior is empty (except for overlap points)
        for h in range(1, height - 1):
            for w in range(1, width - 1):
                if lines[y + h][x + w] not in ' #█':
                    return None
        
        return (x, y, width, height)
    
    rectangles = set()
    # Only start checking from potential top-left corners
    for y in range(height):
        for x in range(width):
            if lines[y][x] in '#█':
                # Check if it's a top-left corner
                if (y == 0 or lines[y-1][x] not in '#█') and (x == 0 or lines[y][x-1] not in '#█'):
                    rect = trace_rectangle(x, y)
                    if rect:
                        rectangles.add(rect)
    
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