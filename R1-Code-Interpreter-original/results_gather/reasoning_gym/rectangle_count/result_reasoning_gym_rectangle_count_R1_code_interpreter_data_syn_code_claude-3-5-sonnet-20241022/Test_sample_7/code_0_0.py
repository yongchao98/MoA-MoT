def find_rectangles(grid):
    # Convert grid to list of strings, removing empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(x, y):
        return grid[y][x] in '#█'
    
    def validate_rectangle(x1, y1, x2, y2):
        # Check if this forms a valid rectangle
        if x1 >= x2 or y1 >= y2:
            return False
            
        # Check all corners
        if not (is_corner(x1, y1) and is_corner(x2, y1) and 
                is_corner(x1, y2) and is_corner(x2, y2)):
            return False
            
        # Check horizontal edges
        for x in range(x1 + 1, x2):
            if not (grid[y1][x] in '#█' and grid[y2][x] in '#█'):
                return False
                
        # Check vertical edges
        for y in range(y1 + 1, y2):
            if not (grid[y][x1] in '#█' and grid[y][x2] in '#█'):
                return False
                
        return True

    rectangles = set()
    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if not is_corner(x1, y1):
                continue
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if validate_rectangle(x1, y1, x2, y2):
                        rectangles.add((x1, y1, x2, y2))

    return len(rectangles)

# The grid as a string (your input grid here)
grid = """                                                                                
                                                                                
                                                                                
                                 ###########################                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
   ##############################█#########                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #########█#################                    
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                            ##########█###################################  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #############################█##########                                  #  
                                ##############################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                            ############################                        
                            #                          #                        
                            ############################                        
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))