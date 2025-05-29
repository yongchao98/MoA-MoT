def find_rectangles(grid):
    lines = [line for line in grid.splitlines() if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def is_complete_rectangle(x1, y1, x2, y2):
        # Check if all corners exist
        if not all(lines[y][x] in '#█' for x, y in [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]):
            return False
            
        # Check horizontal borders
        for x in range(x1 + 1, x2):
            if lines[y1][x] not in '#█' or lines[y2][x] not in '#█':
                return False
                
        # Check vertical borders
        for y in range(y1 + 1, y2):
            if lines[y][x1] not in '#█' or lines[y][x2] not in '#█':
                return False
        
        return True

    rectangles = set()
    # Find top-left corners
    for y1 in range(height-1):
        for x1 in range(width-1):
            if lines[y1][x1] not in '#█':
                continue
                
            # Find matching bottom-right corners
            for y2 in range(y1+1, height):
                for x2 in range(x1+1, width):
                    if lines[y2][x2] not in '#█':
                        continue
                    
                    # Verify if it forms a complete rectangle
                    if is_complete_rectangle(x1, y1, x2, y2):
                        # Store rectangle as tuple of corner coordinates
                        rect = tuple(sorted([(x1, y1), (x2, y2)]))
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