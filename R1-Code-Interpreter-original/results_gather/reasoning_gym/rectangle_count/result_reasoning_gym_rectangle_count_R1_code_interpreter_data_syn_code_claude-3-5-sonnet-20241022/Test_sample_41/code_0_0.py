def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_rectangle(x1, y1, x2, y2):
        # Check if the given coordinates form a valid rectangle
        if x2 <= x1 or y2 <= y1:
            return False
            
        # Check top and bottom borders
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
                
        # Check left and right borders
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
                
        return True

    rectangles = set()
    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] in '#█':
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if grid[y2][x2] in '#█':
                            if is_valid_rectangle(x1, y1, x2, y2):
                                rectangles.add((x1, y1, x2, y2))

    return len(rectangles)

# The ASCII grid
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