def find_rectangles(grid):
    lines = [line for line in grid.split('\n') if line.strip()]
    height = len(lines)
    width = len(lines[0])
    
    def verify_rectangle(top, left, bottom, right):
        # Check if all corners are present
        if not all(lines[y][x] in '#█' for y, x in [
            (top, left), (top, right), 
            (bottom, left), (bottom, right)
        ]):
            return False
            
        # Check horizontal borders
        for x in range(left, right + 1):
            if lines[top][x] not in '#█' or lines[bottom][x] not in '#█':
                return False
                
        # Check vertical borders
        for y in range(top, bottom + 1):
            if lines[y][left] not in '#█' or lines[y][right] not in '#█':
                return False
                
        # Verify it has interior space
        has_space = False
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if lines[y][x] == ' ':
                    has_space = True
                    break
            if has_space:
                break
        
        return has_space

    def find_next_hash(start_y, start_x):
        for y in range(start_y, height):
            x_start = start_x if y == start_y else 0
            for x in range(x_start, width):
                if lines[y][x] in '#█':
                    return (y, x)
        return None

    rectangles = set()
    pos = find_next_hash(0, 0)
    
    while pos:
        top, left = pos
        # Find possible right edges
        for right in range(left + 1, width):
            if lines[top][right] in '#█':
                # Find possible bottom edges
                for bottom in range(top + 1, height):
                    if lines[bottom][right] in '#█' and lines[bottom][left] in '#█':
                        if verify_rectangle(top, left, bottom, right):
                            rectangles.add((top, left, bottom, right))
        
        # Find next starting position
        pos = find_next_hash(top, left + 1)
        if not pos:
            pos = find_next_hash(top + 1, 0)

    return len(rectangles)

grid = """                                                                                
                                                                                
                                                                                
                                #############################                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                         #######█#######################    #                   
                         #      #                      #    #                   
                         #      #                      #    #                   
                         #      #                      #    #                   
                         #######█#######################    #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #############################                   
                                                                                
                                                                                
               ###################################                              
               #                                 #                              
               #                                 #                              
               #                                 #                              
               #                                 #                              
               #                                 #            ###############   
               #                                 #            #             #   
               ###################################            #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                           ###█########     #   
                                                           #  ########█######   
                                                           ############         
                                                                                
                                                                                """

print(find_rectangles(grid))